# Zimbabwe–Zambia Virome Analyses

This repository contains the reproducible R and Python code used for the Zim–Zam microbiome and virome analyses. **No raw data, intermediate files, final results, or figures are included.** Input files must be obtained separately and supplied locally when running the scripts.

## R analyses

The R script is named `viral_profile_characterization.R` and covers the following broad analysis categories:

### Diversity assessment and community composition

Calculation and visualization of within-sample diversity (Shannon diversity) and between-sample community structure using Bray–Curtis dissimilarity, ordination, and group-level comparisons.

### Differential abundance and compositional analyses

Analysis and visualization of differences in microbial, viral, functional, lifestyle, or host-associated feature composition between study groups. The relevant scripts document their filtering, normalization, statistical testing, and multiple-testing procedures.

### Figure generation

Publication figures are generated from locally available analysis results using `ggplot2` and related R packages. Figure scripts do not download or contain study data.

Run an R script from the terminal with:

```bash
Rscript path/to/script.R
```

## Python analyses

### Viral incremental predictive value analysis

The Python script `viral_incremental_prediction.py` evaluates whether viral features improve prediction of fecal metabolite concentrations beyond bacterial features alone.

For each metabolite, the workflow compares:

1. bacterial features only;
2. viral features only; and
3. combined bacterial and viral features.

The analysis uses CLR-transformed microbial abundances, prevalence filtering, training-fold feature selection, Ridge regression, and repeated nested cross-validation. Viral incremental predictive value is calculated as:

```text
Delta R-squared = R-squared(combined model) - R-squared(bacteria-only model)
```

The script accepts the three input files and an output location as command-line arguments. Example:

```bash
python path/to/viral_incremental_prediction.py \
  --three-files-only \
  --tpm-file /path/to/viral_tpm_file.xlsx \
  --bacterial-file /path/to/bacterial_file.txt \
  --metabolite-file /path/to/metabolite_file.xlsx \
  --output-dir /path/to/results
```

The output directory is created locally and is not part of this repository.

## Reproducibility

All file paths should be supplied as command-line arguments or configured locally. Do not hard-code personal directory paths in scripts. Because the repository intentionally contains code only, users must provide authorized input data separately.
