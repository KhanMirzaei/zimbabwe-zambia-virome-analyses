#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"
output_dir <- if (length(args) >= 2) args[[2]] else "results/ancombc2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input <- read_excel(file.path(data_dir, "Zim-Zam.xlsx"), sheet = "df_report")
otu <- input %>%
  select(Contig, starts_with("Zim"), starts_with("Zam")) %>%
  column_to_rownames("Contig") %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ coalesce(.x, 0))) %>%
  as.matrix()
meta <- data.frame(
  Group = ifelse(grepl("^Zim", colnames(otu)), "Zim-R", "Zam-R"),
  row.names = colnames(otu)
)
ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), sample_data(meta))
ps <- prune_taxa(rowSums(otu_table(ps) > 0) >= 2, ps)

fit <- ancombc2(data = ps, fix_formula = "Group", group = "Group",
                p_adj_method = "BH", lib_cut = 0,
                struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05)
res <- as.data.frame(fit$res) %>% rownames_to_column("Contig")
write_tsv(res, file.path(output_dir, "ANCOMBC2_results.tsv"))
