## Author: Mohammadali Khan Mirzaei 

"""Estimate the incremental predictive value of viral features.

from __future__ import annotations

import argparse
import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


METHODS_WORDING = (
    "To assess whether viral features provided predictive information beyond "
    "bacterial composition, metabolite-specific ridge regression models were "
    "fitted using bacterial features alone, viral features alone, or combined "
    "bacterial and viral features. Model performance was evaluated using "
    "repeated nested cross-validation, and the incremental contribution of "
    "viral features was calculated as ΔR² = R²combined − R²bacteria-only. "
    "Mean ΔR² values and 95% confidence intervals were estimated across "
    "repeated cross-validation runs."
)


def log(message: str) -> None:
    print(message, flush=True)


def clean_sample_id(value: object) -> str:
    """Normalize a sample identifier without changing already-clean IDs."""
    label = str(value).strip()
    label = re.sub(r"_R", "", label, flags=re.IGNORECASE)
    label = re.sub(r"^Zi(?!m)", "Zim", label, flags=re.IGNORECASE)
    label = re.sub(r"^Za(?!m)", "Zam", label, flags=re.IGNORECASE)
    return label.strip()


def read_excel_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    frame = pd.read_excel(path)
    if frame.empty:
        raise ValueError(f"Input file is empty: {path}")
    frame = frame.dropna(axis=0, how="all").dropna(axis=1, how="all")
    frame.columns = frame.columns.map(lambda x: str(x).strip())
    frame = frame.loc[:, frame.columns != ""]
    return frame


def read_otu_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    frame = pd.read_csv(path, sep="\t")
    if frame.empty or frame.shape[1] < 2:
        raise ValueError(f"OTU input is empty or malformed: {path}")
    frame = frame.dropna(axis=0, how="all").dropna(axis=1, how="all")
    frame.columns = frame.columns.map(lambda x: str(x).strip())
    return frame


def numeric_and_deduplicated(
    frame: pd.DataFrame, *, microbial: bool
) -> tuple[pd.DataFrame, int, int]:
    """Clean a sample x feature table and resolve duplicate labels."""
    result = frame.copy()
    result.index = result.index.map(clean_sample_id)
    result.columns = result.columns.map(lambda x: str(x).strip())
    duplicate_samples = int(result.index.duplicated(keep=False).sum())
    duplicate_features = int(result.columns.duplicated(keep=False).sum())
    result = result.apply(pd.to_numeric, errors="coerce")
    result = result.replace([np.inf, -np.inf], np.nan)

    if result.index.has_duplicates:
        result = result.groupby(level=0, sort=False).mean()
    if result.columns.has_duplicates:
        transposed = result.T
        if microbial:
            result = transposed.groupby(level=0, sort=False).sum(min_count=1).T
        else:
            result = transposed.groupby(level=0, sort=False).mean().T

    result = result.dropna(axis=1, how="all")
    if microbial:
        negative = (result < 0).any(axis=0)
        if negative.any():
            examples = ", ".join(map(str, result.columns[negative][:5]))
            raise ValueError(
                "CLR requires non-negative microbial abundances; negative values "
                f"were found in: {examples}"
            )
    return result, duplicate_samples, duplicate_features


def prepare_inputs(
    virome_path: Path,
    tpm_path: Path,
    metabolite_path: Path,
    otu_path: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Validate raw uploads, match contigs/samples, and save intermediates."""
    log("Reading and validating the four uploaded source files ...")
    virome_raw = read_excel_table(virome_path)
    tpm_raw = read_excel_table(tpm_path)
    metabolite_raw = read_excel_table(metabolite_path)
    otu_raw = read_otu_table(otu_path)

    qc: list[dict[str, object]] = []

    def record(metric: str, value: object, status: str, detail: str = "") -> None:
        qc.append({"Metric": metric, "Value": value, "Status": status, "Detail": detail})

    record("virome_rows", len(virome_raw), "PASS")
    record("tpm_rows", len(tpm_raw), "PASS")
    record("metabolite_samples_raw", len(metabolite_raw), "PASS")
    record("bacterial_features_raw", len(otu_raw), "PASS")

    virome_id = virome_raw.columns[0]
    tpm_id = tpm_raw.columns[0]
    metabolite_id = metabolite_raw.columns[0]
    otu_id = otu_raw.columns[0]
    virome_raw[virome_id] = virome_raw[virome_id].map(
        lambda x: "" if pd.isna(x) else str(x).strip()
    )
    tpm_raw[tpm_id] = tpm_raw[tpm_id].map(
        lambda x: "" if pd.isna(x) else str(x).strip()
    )
    virome_raw = virome_raw.loc[virome_raw[virome_id] != ""].copy()
    tpm_raw = tpm_raw.loc[tpm_raw[tpm_id] != ""].copy()

    virome_dup = int(virome_raw[virome_id].duplicated(keep=False).sum())
    tpm_dup = int(tpm_raw[tpm_id].duplicated(keep=False).sum())
    record(
        "virome_duplicate_contig_occurrences",
        virome_dup,
        "PASS" if virome_dup == 0 else "WARN",
        "First annotation row is retained for duplicated contig IDs.",
    )
    record(
        "tpm_duplicate_contig_occurrences",
        tpm_dup,
        "PASS" if tpm_dup == 0 else "WARN",
        "Duplicated abundance rows are summed after sample matching.",
    )
    virome_unique = virome_raw.drop_duplicates(virome_id, keep="first").set_index(virome_id)
    main_contigs = set(virome_unique.index)
    tpm_contigs = set(tpm_raw[tpm_id])
    shared_contigs = tpm_contigs & main_contigs
    record("tpm_unique_contigs", len(tpm_contigs), "PASS")
    record("main_virome_unique_contigs", len(main_contigs), "PASS")
    record(
        "tpm_contigs_in_main_virome",
        len(shared_contigs),
        "PASS" if len(shared_contigs) == len(tpm_contigs) else "WARN",
    )
    record(
        "tpm_contigs_excluded_not_in_main_virome",
        len(tpm_contigs - main_contigs),
        "PASS" if not (tpm_contigs - main_contigs) else "WARN",
        "Only contigs present in both files are retained.",
    )
    if not shared_contigs:
        raise ValueError("No TPM contigs match the authoritative virome table.")
    tpm_raw = tpm_raw.loc[tpm_raw[tpm_id].isin(shared_contigs)].copy()

    tpm_values = tpm_raw.drop(columns=tpm_id).apply(pd.to_numeric, errors="coerce")
    tpm_missing = int(tpm_values.isna().sum().sum())
    tpm_negative = int((tpm_values < 0).sum().sum())
    tpm_all_zero_features = int((tpm_values.fillna(0).sum(axis=1) == 0).sum())
    tpm_all_zero_samples = int((tpm_values.fillna(0).sum(axis=0) == 0).sum())
    tpm_totals = tpm_values.sum(axis=0)
    tpm_is_million = bool(np.allclose(tpm_totals, 1_000_000, rtol=1e-6, atol=1.0))
    detected_per_sample = (tpm_values.fillna(0) > 0).sum(axis=0)
    max_feature_fraction = tpm_values.max(axis=0).div(tpm_totals.replace(0, np.nan))
    highly_dominated = max_feature_fraction[max_feature_fraction > 0.90]
    record("tpm_missing_or_nonnumeric_values", tpm_missing, "PASS" if tpm_missing == 0 else "WARN")
    record("tpm_negative_values", tpm_negative, "PASS" if tpm_negative == 0 else "FAIL")
    record("tpm_all_zero_contigs", tpm_all_zero_features, "PASS" if tpm_all_zero_features == 0 else "WARN")
    record("tpm_all_zero_samples", tpm_all_zero_samples, "PASS" if tpm_all_zero_samples == 0 else "FAIL")
    record(
        "tpm_sample_sums_equal_1e6",
        tpm_is_million,
        "PASS" if tpm_is_million else "WARN",
        f"Range: {tpm_totals.min():.6g} to {tpm_totals.max():.6g}",
    )
    record(
        "tpm_detected_contigs_per_sample_range",
        f"{int(detected_per_sample.min())}-{int(detected_per_sample.max())}",
        "INFO",
    )
    record(
        "tpm_max_single_contig_fraction_range",
        f"{max_feature_fraction.min():.6g}-{max_feature_fraction.max():.6g}",
        "WARN" if len(highly_dominated) else "INFO",
        "Samples above 90%: " + ", ".join(highly_dominated.index.astype(str)),
    )
    if tpm_negative or tpm_all_zero_samples:
        raise ValueError("TPM validation failed; see ml_input_qc.csv for details.")

    virus_oriented = tpm_raw.set_index(tpm_id).T
    virus, v_dup_s, v_dup_f = numeric_and_deduplicated(virus_oriented, microbial=True)

    ignored_otu_columns = {
        col
        for col in otu_raw.columns[1:]
        if str(col).strip().lower() == "taxonomy"
        or str(col).strip().lower().startswith("unnamed")
    }
    bacteria_oriented = otu_raw.drop(columns=list(ignored_otu_columns)).set_index(otu_id).T
    bacteria, b_dup_s, b_dup_f = numeric_and_deduplicated(
        bacteria_oriented, microbial=True
    )

    metabolites_oriented = metabolite_raw.set_index(metabolite_id)
    metabolites, m_dup_s, m_dup_f = numeric_and_deduplicated(
        metabolites_oriented, microbial=False
    )
    record("viral_duplicate_cleaned_samples", v_dup_s, "PASS" if v_dup_s == 0 else "WARN")
    record("bacterial_duplicate_cleaned_samples", b_dup_s, "PASS" if b_dup_s == 0 else "WARN")
    record("metabolite_duplicate_cleaned_samples", m_dup_s, "PASS" if m_dup_s == 0 else "WARN")
    record("viral_duplicate_features", v_dup_f, "PASS" if v_dup_f == 0 else "WARN")
    record("bacterial_duplicate_features", b_dup_f, "PASS" if b_dup_f == 0 else "WARN")
    record("duplicate_metabolites", m_dup_f, "PASS" if m_dup_f == 0 else "WARN")

    all_samples = sorted(set(bacteria.index) | set(virus.index) | set(metabolites.index))
    shared_samples = sorted(set(bacteria.index) & set(virus.index) & set(metabolites.index))
    if len(shared_samples) < 4:
        raise ValueError(f"Only {len(shared_samples)} shared samples; at least 4 are required.")
    sample_report = pd.DataFrame(
        {
            "SampleID": all_samples,
            "In_bacteria": [sample in bacteria.index for sample in all_samples],
            "In_viral_TPM": [sample in virus.index for sample in all_samples],
            "In_metabolites": [sample in metabolites.index for sample in all_samples],
            "Included_in_ML": [sample in shared_samples for sample in all_samples],
        }
    )
    bacteria = bacteria.loc[shared_samples]
    virus = virus.loc[shared_samples]
    metabolites = metabolites.loc[shared_samples]
    record("shared_samples_for_ml", len(shared_samples), "PASS")
    record("bacterial_features_matched", bacteria.shape[1], "PASS")
    record("viral_features_matched", virus.shape[1], "PASS")
    record("metabolites_matched", metabolites.shape[1], "PASS")

    required = max(1, min(3, int(np.ceil(0.10 * len(shared_samples)))))
    record("full_data_prevalence_threshold", required, "INFO", "Model filtering is refit inside each CV training split.")
    record("bacterial_features_passing_full_data_prevalence", int(((bacteria.fillna(0) > 0).sum(axis=0) >= required).sum()), "INFO")
    record("viral_features_passing_full_data_prevalence", int(((virus.fillna(0) > 0).sum(axis=0) >= required).sum()), "INFO")

    # Keep only non-abundance metadata from the main virome workbook.
    viral_sample_labels = {clean_sample_id(col) for col in tpm_raw.columns if col != tpm_id}
    abundance_columns = [
        col for col in virome_unique.columns if clean_sample_id(col) in viral_sample_labels
    ]
    annotation_columns = [col for col in virome_unique.columns if col not in abundance_columns]
    annotations = virome_unique.reindex(virus.columns)[annotation_columns]

    output_dir.mkdir(parents=True, exist_ok=True)
    bacteria.to_csv(output_dir / "ml_intermediate_bacteria_matched.csv", index_label="SampleID")
    virus.to_csv(output_dir / "ml_intermediate_viral_tpm_matched.csv", index_label="SampleID")
    metabolites.to_csv(output_dir / "ml_intermediate_metabolites_matched.csv", index_label="SampleID")
    annotations.to_csv(output_dir / "ml_intermediate_virome_annotations_matched.csv", index_label="Contig")
    sample_report.to_csv(output_dir / "ml_sample_matching_report.csv", index=False)
    qc_df = pd.DataFrame(qc)
    qc_df.to_csv(output_dir / "ml_input_qc.csv", index=False)

    log(
        f"Matched {len(shared_samples)} samples, {bacteria.shape[1]} bacterial "
        f"features, {virus.shape[1]} viral features, and {metabolites.shape[1]} metabolites."
    )
    log("Intermediate matched files and QC reports were written.")
    return bacteria, virus, metabolites, qc_df, sample_report


def prepare_three_inputs(
    tpm_path: Path,
    metabolite_path: Path,
    otu_path: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load, orient, match, and save intermediates using three files only."""
    log("Reading the three uploaded files (no virome annotation workbook required) ...")
    tpm_raw = read_excel_table(tpm_path)
    metabolite_raw = read_excel_table(metabolite_path)
    otu_raw = read_otu_table(otu_path)

    tpm_id = tpm_raw.columns[0]
    tpm = numeric_and_deduplicated(tpm_raw.set_index(tpm_id).T, microbial=True)[0]
    ignored = {
        col for col in otu_raw.columns[1:]
        if str(col).strip().lower() == "taxonomy"
        or str(col).strip().lower().startswith("unnamed")
    }
    otu_id = otu_raw.columns[0]
    bacteria = numeric_and_deduplicated(
        otu_raw.drop(columns=list(ignored)).set_index(otu_id).T,
        microbial=True,
    )[0]

    metabolite_id = metabolite_raw.columns[0]
    m = metabolite_raw.set_index(metabolite_id)
    m_feature_by_sample = m.T
    m_sample_by_feature = m
    microbial_samples = set(tpm.index) & set(bacteria.index)
    candidate_a = set(map(clean_sample_id, m_feature_by_sample.index)) & microbial_samples
    candidate_b = set(map(clean_sample_id, m_sample_by_feature.index)) & microbial_samples
    m_oriented = m_feature_by_sample if len(candidate_a) >= len(candidate_b) else m_sample_by_feature
    metabolites = numeric_and_deduplicated(m_oriented, microbial=False)[0]

    shared = sorted(set(tpm.index) & set(bacteria.index) & set(metabolites.index))
    if len(shared) < 4:
        raise ValueError(f"Only {len(shared)} shared samples; at least 4 are required.")
    tpm, bacteria, metabolites = tpm.loc[shared], bacteria.loc[shared], metabolites.loc[shared]
    output_dir.mkdir(parents=True, exist_ok=True)
    tpm.to_csv(output_dir / "ml_intermediate_viral_tpm_matched.csv", index_label="SampleID")
    bacteria.to_csv(output_dir / "ml_intermediate_bacteria_matched.csv", index_label="SampleID")
    metabolites.to_csv(output_dir / "ml_intermediate_metabolites_matched.csv", index_label="SampleID")
    pd.DataFrame({"SampleID": shared, "Included_in_ML": True}).to_csv(
        output_dir / "ml_sample_matching_report.csv", index=False
    )
    qc = pd.DataFrame([
        {"Metric": "shared_samples", "Value": len(shared), "Status": "PASS"},
        {"Metric": "viral_features", "Value": tpm.shape[1], "Status": "PASS"},
        {"Metric": "bacterial_features", "Value": bacteria.shape[1], "Status": "PASS"},
        {"Metric": "metabolites", "Value": metabolites.shape[1], "Status": "PASS"},
        {"Metric": "tpm_negative_values", "Value": int((tpm < 0).sum().sum()), "Status": "PASS"},
    ])
    qc.to_csv(output_dir / "ml_input_qc.csv", index=False)
    log(f"Matched {len(shared)} samples, {bacteria.shape[1]} bacterial features, "
        f"{tpm.shape[1]} viral features, and {metabolites.shape[1]} metabolites.")
    return bacteria, tpm, metabolites, qc, pd.DataFrame({"SampleID": shared})


class BlockwisePrevalenceCLR(BaseEstimator, TransformerMixin):
    """Training-only prevalence filter and separate CLR per feature block."""

    def __init__(
        self,
        block_sizes: tuple[int, ...],
        min_samples: int = 3,
        min_fraction: float = 0.10,
        max_features: int | None = None,
    ) -> None:
        self.block_sizes = block_sizes
        self.min_samples = min_samples
        self.min_fraction = min_fraction
        self.max_features = max_features

    def fit(self, X: np.ndarray, y: np.ndarray | None = None):
        values = np.asarray(X, dtype=float)
        if values.ndim != 2 or values.shape[1] != sum(self.block_sizes):
            raise ValueError("Unexpected feature matrix shape in CLR transformer.")
        if np.nanmin(values) < 0:
            raise ValueError("CLR input contains negative values.")

        # "At least 3 samples OR at least 10%" is the less restrictive count.
        required = max(
            1,
            min(int(self.min_samples), int(np.ceil(self.min_fraction * len(values)))),
        )
        self.required_presence_ = required
        self.masks_: list[np.ndarray] = []
        self.pseudocounts_: list[float] = []
        start = 0
        for size in self.block_sizes:
            block = values[:, start : start + size]
            finite = np.where(np.isfinite(block), block, 0.0)
            mask = np.sum(finite > 0, axis=0) >= required
            if not np.any(mask):
                # Preserve one column so an extremely sparse block degrades to
                # an intercept-only contribution instead of crashing.
                mask[int(np.argmax(np.sum(finite > 0, axis=0)))] = True
            kept = finite[:, mask]
            positive = kept[kept > 0]
            pseudocount = float(0.5 * positive.min()) if positive.size else 1.0
            if self.max_features is not None and kept.shape[1] > self.max_features:
                logged = np.log(np.where(kept > 0, kept, pseudocount))
                clr = logged - logged.mean(axis=1, keepdims=True)
                variances = np.nanvar(clr, axis=0)
                keep_local = np.argsort(variances)[-int(self.max_features):]
                reduced = np.zeros_like(mask)
                original_indices = np.flatnonzero(mask)
                reduced[original_indices[keep_local]] = True
                mask = reduced
            self.masks_.append(mask)
            self.pseudocounts_.append(pseudocount)
            start += size
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        values = np.asarray(X, dtype=float)
        transformed: list[np.ndarray] = []
        start = 0
        for size, mask, pseudocount in zip(
            self.block_sizes, self.masks_, self.pseudocounts_
        ):
            block = values[:, start : start + size][:, mask]
            block = np.where(np.isfinite(block), block, 0.0)
            block = np.where(block > 0, block, pseudocount)
            logged = np.log(block)
            transformed.append(logged - logged.mean(axis=1, keepdims=True))
            start += size
        return np.concatenate(transformed, axis=1)


def make_pipeline(
    block_sizes: tuple[int, ...], max_features_per_block: int | None = None
) -> Pipeline:
    return Pipeline(
        steps=[
            ("clr", BlockwisePrevalenceCLR(
                block_sizes=block_sizes,
                max_features=max_features_per_block,
            )),
            ("scale", StandardScaler()),
            ("ridge", Ridge(solver="lsqr", max_iter=20_000)),
        ]
    )


def fit_predict_nested(
    X: np.ndarray,
    y: np.ndarray,
    train: np.ndarray,
    test: np.ndarray,
    block_sizes: tuple[int, ...],
    inner_splits: list[tuple[np.ndarray, np.ndarray]],
    alphas: np.ndarray,
    max_features_per_block: int | None,
) -> tuple[np.ndarray, float]:
    search = GridSearchCV(
        estimator=make_pipeline(block_sizes, max_features_per_block),
        param_grid={"ridge__alpha": alphas},
        scoring="neg_mean_squared_error",
        cv=inner_splits,
        n_jobs=1,
        refit=True,
        error_score="raise",
    )
    search.fit(X[train], y[train])
    prediction = search.predict(X[test])
    return prediction, float(search.best_params_["ridge__alpha"])


def bootstrap_mean_ci(
    values: np.ndarray, n_bootstrap: int, rng: np.random.Generator
) -> tuple[float, float]:
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if clean.size == 0:
        return np.nan, np.nan
    if clean.size == 1:
        return float(clean[0]), float(clean[0])
    # Chunking avoids a large allocation if users request many repeats/bootstraps.
    boot_means = np.empty(n_bootstrap, dtype=float)
    chunk = 2_000
    for start in range(0, n_bootstrap, chunk):
        stop = min(start + chunk, n_bootstrap)
        draws = rng.choice(clean, size=(stop - start, clean.size), replace=True)
        boot_means[start:stop] = draws.mean(axis=1)
    low, high = np.quantile(boot_means, [0.025, 0.975])
    return float(low), float(high)


def analyze_metabolite(
    name: str,
    y_series: pd.Series,
    bacteria: pd.DataFrame,
    virus: pd.DataFrame,
    repeats: int,
    outer_folds_requested: int,
    inner_folds_requested: int,
    alphas: np.ndarray,
    seed: int,
    max_features_per_block: int | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    valid = y_series.notna()
    y = y_series.loc[valid].to_numpy(dtype=float)
    b = bacteria.loc[valid].to_numpy(dtype=float)
    v = virus.loc[valid].to_numpy(dtype=float)
    n = len(y)

    base_summary: dict[str, object] = {"Metabolite": name, "N_samples": n}
    if n < 4:
        return [], {
            **base_summary,
            "N_repeats": 0,
            "Mean_R2_bacteria_only": np.nan,
            "Mean_R2_virus_only": np.nan,
            "Mean_R2_combined": np.nan,
            "Mean_Delta_R2": np.nan,
            "Median_Delta_R2": np.nan,
            "SD_Delta_R2": np.nan,
            "CI95_lower": np.nan,
            "CI95_upper": np.nan,
            "Classification": "Uncertain",
            "Status": "skipped_too_few_samples",
        }
    if np.nanvar(y) <= np.finfo(float).eps:
        return [], {
            **base_summary,
            "N_repeats": 0,
            "Mean_R2_bacteria_only": np.nan,
            "Mean_R2_virus_only": np.nan,
            "Mean_R2_combined": np.nan,
            "Mean_Delta_R2": np.nan,
            "Median_Delta_R2": np.nan,
            "SD_Delta_R2": np.nan,
            "CI95_lower": np.nan,
            "CI95_upper": np.nan,
            "Classification": "Uncertain",
            "Status": "skipped_constant_metabolite",
        }

    outer_folds = min(outer_folds_requested, n)
    rows: list[dict[str, object]] = []
    combined = np.concatenate([b, v], axis=1)
    model_specs = {
        "bacteria": (b, (b.shape[1],)),
        "virus": (v, (v.shape[1],)),
        "combined": (combined, (b.shape[1], v.shape[1])),
    }

    for repeat in range(repeats):
        outer = KFold(n_splits=outer_folds, shuffle=True, random_state=seed + repeat)
        outer_split_list = list(outer.split(np.arange(n)))
        predictions = {key: np.full(n, np.nan) for key in model_specs}
        selected_alphas = {key: [] for key in model_specs}
        inner_counts: list[int] = []

        for fold, (train, test) in enumerate(outer_split_list):
            inner_folds = min(inner_folds_requested, len(train))
            if inner_folds < 2:
                raise ValueError("Too few outer-training samples for inner CV.")
            inner = KFold(
                n_splits=inner_folds,
                shuffle=True,
                random_state=seed + repeat * 1_000 + fold,
            )
            # Materializing once guarantees identical inner splits for all models.
            inner_splits = list(inner.split(np.arange(len(train))))
            inner_counts.append(inner_folds)
            for key, (X, block_sizes) in model_specs.items():
                pred, alpha = fit_predict_nested(
                    X, y, train, test, block_sizes, inner_splits, alphas,
                    max_features_per_block,
                )
                predictions[key][test] = pred
                selected_alphas[key].append(alpha)

        scores = {key: r2_score(y, pred) for key, pred in predictions.items()}
        delta = scores["combined"] - scores["bacteria"]
        rows.append(
            {
                "Metabolite": name,
                "Repeat": repeat + 1,
                "N_samples": n,
                "Outer_folds": outer_folds,
                "Inner_folds_min": min(inner_counts),
                "R2_bacteria_only": scores["bacteria"],
                "R2_virus_only": scores["virus"],
                "R2_combined": scores["combined"],
                "Delta_R2": delta,
                "Mean_selected_alpha_bacteria": np.mean(selected_alphas["bacteria"]),
                "Mean_selected_alpha_virus": np.mean(selected_alphas["virus"]),
                "Mean_selected_alpha_combined": np.mean(selected_alphas["combined"]),
                "Status": "ok",
            }
        )

    result = pd.DataFrame(rows)
    delta_values = result["Delta_R2"].to_numpy()
    rng = np.random.default_rng(seed + 100_000)
    low, high = bootstrap_mean_ci(delta_values, 10_000, rng)
    if low > 0:
        classification = "Added predictive value"
    elif high < 0:
        classification = "Reduced predictive value"
    else:
        classification = "Uncertain"

    summary = {
        **base_summary,
        "N_repeats": repeats,
        "Mean_R2_bacteria_only": result["R2_bacteria_only"].mean(),
        "Mean_R2_virus_only": result["R2_virus_only"].mean(),
        "Mean_R2_combined": result["R2_combined"].mean(),
        "Mean_Delta_R2": result["Delta_R2"].mean(),
        "Median_Delta_R2": result["Delta_R2"].median(),
        "SD_Delta_R2": result["Delta_R2"].std(ddof=1),
        "CI95_lower": low,
        "CI95_upper": high,
        "Classification": classification,
        "Status": "ok",
    }
    return rows, summary


def display_label(name: str) -> str:
    return name.replace("_", " ").replace("µmol/g", "(µmol/g)")


def create_plot(summary: pd.DataFrame, png_path: Path, pdf_path: Path) -> None:
    plot_data = summary.loc[summary["Status"].eq("ok")].copy()
    plot_data = plot_data.sort_values("Mean_Delta_R2", ascending=False)
    if plot_data.empty:
        raise ValueError("No analyzable metabolites were available for plotting.")

    palette = {
        "Added predictive value": "#2E8B57",
        "Reduced predictive value": "#C44E52",
        "Uncertain": "#7A7A7A",
    }
    x = np.arange(len(plot_data))
    means = plot_data["Mean_Delta_R2"].to_numpy()
    lower = plot_data["CI95_lower"].to_numpy()
    upper = plot_data["CI95_upper"].to_numpy()
    yerr = np.vstack([means - lower, upper - means])
    colors = plot_data["Classification"].map(palette).to_list()

    width = max(11.0, 0.48 * len(plot_data))
    fig, ax = plt.subplots(figsize=(width, 6.7))
    ax.bar(
        x,
        means,
        yerr=yerr,
        color=colors,
        edgecolor="black",
        linewidth=0.45,
        error_kw={"ecolor": "black", "elinewidth": 1.0, "capsize": 3},
    )
    ax.axhline(0, color="black", linewidth=1.0)
    ax.set_ylabel("Mean ΔR² (combined − bacteria only)", fontsize=11)
    ax.set_xlabel("Metabolite", fontsize=11)
    ax.set_title("Incremental predictive value of viral features", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [display_label(value) for value in plot_data["Metabolite"]],
        rotation=60,
        ha="right",
        fontsize=8.5,
    )
    ax.grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    present = [label for label in palette if label in set(plot_data["Classification"])]
    handles = [Patch(facecolor=palette[label], edgecolor="black", label=label) for label in present]
    ax.legend(handles=handles, frameon=False, loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--virome-file", default="Zim-Zam.xlsx")
    parser.add_argument("--tpm-file", default="TPM_filtered_nonzero_contigs.xlsx")
    parser.add_argument("--metabolite-file", default="Zim_Zam_SCFA_BA_v2.xlsx")
    parser.add_argument("--otu-file", default="OTUs-Table_Zim_Zam.txt")
    parser.add_argument(
        "--three-files-only", action="store_true",
        help="Use TPM, bacterial, and metabolite files without Zim-Zam.xlsx.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Validate and write matched intermediates without fitting models.",
    )
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--alpha-min-log10", type=float, default=-4.0)
    parser.add_argument("--alpha-max-log10", type=float, default=4.0)
    parser.add_argument("--n-alphas", type=int, default=17)
    parser.add_argument(
        "--max-features-per-block",
        type=int,
        default=100,
        help="Training-fold-only variance reduction per bacterial/viral block (default: 100).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.outer_folds < 2 or args.inner_folds < 2 or args.repeats < 1:
        raise ValueError("CV folds must be at least 2 and repeats at least 1.")
    if args.bootstrap < 1 or args.n_alphas < 2:
        raise ValueError("Bootstrap count must be positive and n-alphas at least 2.")

    data_dir = args.data_dir.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve() if args.output_dir is not None else data_dir
    )
    log(f"Data directory: {data_dir}")
    log(f"Output directory: {output_dir}")
    if args.max_features_per_block is None:
        log("Feature preselection: OFF (all prevalence-passing features)")
    else:
        log(
            "Feature preselection: ON (top "
            f"{args.max_features_per_block} variance-ranked features per block, "
            "selected inside each training split)"
        )
    if args.three_files_only:
        bacteria, virus, metabolites, _, _ = prepare_three_inputs(
            data_dir / args.tpm_file,
            data_dir / args.metabolite_file,
            data_dir / args.otu_file,
            output_dir,
        )
    else:
        bacteria, virus, metabolites, _, _ = prepare_inputs(
            data_dir / args.virome_file,
            data_dir / args.tpm_file,
            data_dir / args.metabolite_file,
            data_dir / args.otu_file,
            output_dir,
        )
    log(
        "Missing microbial entries will be treated as zeros inside the CLR "
        f"transformer: bacteria={int(bacteria.isna().sum().sum())}, "
        f"virus={int(virus.isna().sum().sum())}."
    )
    if args.prepare_only:
        log("Preparation-only run completed; model fitting was not requested.")
        return 0

    alphas = np.logspace(args.alpha_min_log10, args.alpha_max_log10, args.n_alphas)
    all_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    log(
        f"Starting nested CV: {args.outer_folds} outer folds x {args.repeats} "
        f"repeats; up to {args.inner_folds} inner folds; {len(alphas)} alphas."
    )
    for number, metabolite in enumerate(metabolites.columns, start=1):
        log(f"[{number:02d}/{metabolites.shape[1]:02d}] {metabolite}")
        rows, summary = analyze_metabolite(
            metabolite,
            metabolites[metabolite],
            bacteria,
            virus,
            repeats=args.repeats,
            outer_folds_requested=args.outer_folds,
            inner_folds_requested=args.inner_folds,
            alphas=alphas,
            seed=args.seed + number * 10_000,
            max_features_per_block=args.max_features_per_block,
        )
        # Recompute the requested bootstrap size here, allowing a CLI override.
        if rows:
            deltas = np.array([float(row["Delta_R2"]) for row in rows])
            low, high = bootstrap_mean_ci(
                deltas,
                args.bootstrap,
                np.random.default_rng(args.seed + number * 100_000),
            )
            summary["CI95_lower"] = low
            summary["CI95_upper"] = high
            summary["Classification"] = (
                "Added predictive value"
                if low > 0
                else "Reduced predictive value"
                if high < 0
                else "Uncertain"
            )
            log(
                f"    mean ΔR²={summary['Mean_Delta_R2']:.4f}; "
                f"95% CI [{low:.4f}, {high:.4f}]; {summary['Classification']}"
            )
        else:
            log(f"    {summary['Status']}")
        all_rows.extend(rows)
        summaries.append(summary)

    repeats_df = pd.DataFrame(all_rows)
    summary_df = pd.DataFrame(summaries).sort_values(
        "Mean_Delta_R2", ascending=False, na_position="last"
    )
    repeat_path = output_dir / "ml_phage_incremental_r2_repeats.csv"
    summary_path = output_dir / "ml_phage_incremental_r2_summary.csv"
    png_path = output_dir / "ml_phage_incremental_r2_plot.png"
    pdf_path = output_dir / "ml_phage_incremental_r2_plot.pdf"
    repeats_df.to_csv(repeat_path, index=False, float_format="%.10g")
    summary_df.to_csv(summary_path, index=False, float_format="%.10g")
    create_plot(summary_df, png_path, pdf_path)

    log("\nCompleted successfully. Output files:")
    for path in (repeat_path, summary_path, png_path, pdf_path):
        log(f"  {path}")
    log("\nMethods wording:\n" + METHODS_WORDING)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        warnings.warn(f"Analysis failed: {exc}")
        raise
