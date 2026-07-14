#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"
output_dir <- if (length(args) >= 2) args[[2]] else "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tpm <- read_excel(file.path(data_dir, "TPM_filtered_nonzero_contigs.xlsx"))
annotation <- read_excel(file.path(data_dir, "filtered_nonzero_contigs.xlsx"))
sample_cols <- grep("^(Zim|Zam)", names(tpm), value = TRUE)

long <- tpm %>%
  select(Contig, all_of(sample_cols)) %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "TPM") %>%
  mutate(TPM = replace_na(as.numeric(TPM), 0),
         Group = ifelse(grepl("^Zim", Sample), "Zim-R", "Zam-R")) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = TPM / sum(TPM)) %>%
  ungroup() %>%
  left_join(annotation, by = "Contig")

plot_composition <- function(variable, title, file_name) {
  dat <- long %>%
    filter(!is.na(.data[[variable]]), .data[[variable]] != "") %>%
    group_by(Sample, Group, Category = .data[[variable]]) %>%
    summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop")
  p <- ggplot(dat, aes(Category, RelativeAbundance, fill = Group)) +
    geom_boxplot(position = position_dodge(0.7), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.12), size = 1.5) +
    scale_fill_manual(values = c("Zim-R" = "#4E79A7", "Zam-R" = "#59A14F")) +
    labs(title = title, x = NULL, y = "Relative abundance", fill = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_dir, file_name), p, width = 7, height = 5, dpi = 300)
}

plot_composition("lifestyle", "Replication-cycle composition",
                 "replication_cycle_composition.png")
plot_composition("vhost_hostPhylum", "Predicted host-phylum composition",
                 "predicted_host_phylum_composition.png")
