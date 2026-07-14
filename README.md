# Zimbabwe–Zambia Virome Analyses

This repository contains the reproducible R and Python code used for the virome analyses in the Zimbabwe–Zambia manuscript. It contains **code only**. No raw data, intermediate files, result tables, or figures are included.

## Repository structure

```text
R/
  01_virome_diversity.R
  02_replication_cycle_host_characterization.R
  03_virome_differential_abundance_ANCOMBC2.R
Python/
  viral_incremental_prediction.py
```

## Local input files

Obtain the authorized input files separately and place them in a local data directory. Expected filenames are:

```text
TPM_filtered_nonzero_contigs.xlsx
filtered_nonzero_contigs.xlsx
Zim-Zam.xlsx
OTUs-Table_Zim_Zam.txt
Zim_Zam_SCFA_BA_v2.xlsx
```

The scripts never require the data to be copied into this repository.

## R analyses

### Virome diversity and community composition

`R/01_virome_diversity.R` calculates Shannon diversity, Bray–Curtis dissimilarity, PERMANOVA, and a two-dimensional t-SNE visualization.

Required input: `TPM_filtered_nonzero_contigs.xlsx`.

Run:

```bash
Rscript R/01_virome_diversity.R /path/to/data /path/to/results
```

### Replication-cycle and predicted-host characterization

`R/02_replication_cycle_host_characterization.R` summarizes viral relative abundance by replication-cycle category and predicted host phylum.

Required inputs: `TPM_filtered_nonzero_contigs.xlsx` and `filtered_nonzero_contigs.xlsx`.

Run:

```bash
Rscript R/02_replication_cycle_host_characterization.R /path/to/data /path/to/results
```

### Virome differential-abundance analysis

`R/03_virome_differential_abundance_ANCOMBC2.R` performs group-level viral differential-abundance testing with ANCOM-BC2 and writes a tab-delimited result table.

Required input: `Zim-Zam.xlsx`, containing the `df_report` sheet.

Run:

```bash
Rscript R/03_virome_differential_abundance_ANCOMBC2.R /path/to/data /path/to/results
```

R packages used include `readxl`, `dplyr`, `tidyr`, `vegan`, `Rtsne`, `ggplot2`, `phyloseq`, and `ANCOMBC`.

## Python analysis

### Viral incremental predictive value

`Python/viral_incremental_prediction.py` evaluates whether viral features improve prediction of fecal metabolite concentrations beyond bacterial features alone. It compares bacterial-only, viral-only, and combined Ridge regression models using CLR-transformed microbial abundances, prevalence filtering, training-fold feature selection, and repeated nested cross-validation.

The incremental metric is:

```text
Delta R-squared = R-squared(combined model) - R-squared(bacteria-only model)
```

Required inputs: `TPM_filtered_nonzero_contigs.xlsx`, `OTUs-Table_Zim_Zam.txt`, and `Zim_Zam_SCFA_BA_v2.xlsx`.

Run:

```bash
python Python/viral_incremental_prediction.py \
  --three-files-only \
  --data-dir /path/to/data \
  --output-dir /path/to/results
```

Python packages used include `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `openpyxl`, and `scipy`.

## Reproducibility and data availability

All paths are supplied locally at run time; no personal computer paths are embedded in the workflow. The underlying study data are not distributed with this repository and must be accessed through the appropriate data-sharing procedure.
