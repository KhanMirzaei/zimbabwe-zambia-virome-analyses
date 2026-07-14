# Zimbabwe–Zambia Virome Analyses

Code used for the virome analyses in the Zimbabwe–Zambia manuscript. No data, intermediate files, results, or figures are included.

## Files

```text
R/01_virome_diversity.R
R/02_replication_cycle_host_characterization.R
R/03_virome_differential_abundance_ANCOMBC2.R
Python/viral_incremental_prediction.py
```

## Input files

Place authorized input files in a local data directory:

```text
TPM_filtered_nonzero_contigs.xlsx
filtered_nonzero_contigs.xlsx
Zim-Zam.xlsx
OTUs-Table_Zim_Zam.txt
Zim_Zam_SCFA_BA_v2.xlsx
```

## Run the analyses

```bash
Rscript R/01_virome_diversity.R /path/to/data /path/to/results
Rscript R/02_replication_cycle_host_characterization.R /path/to/data /path/to/results
Rscript R/03_virome_differential_abundance_ANCOMBC2.R /path/to/data /path/to/results
```

```bash
python Python/viral_incremental_prediction.py \
  --three-files-only \
  --data-dir /path/to/data \
  --output-dir /path/to/results
```

R analyses use `vegan`, `Rtsne`, `ggplot2`, `phyloseq`, and `ANCOMBC`. The Python workflow uses `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `openpyxl`, and `scipy`.
