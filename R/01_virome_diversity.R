#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(vegan)
  library(ggplot2)
  library(Rtsne)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"
output_dir <- if (length(args) >= 2) args[[2]] else "results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tpm <- read_excel(file.path(data_dir, "TPM_filtered_nonzero_contigs.xlsx"))
sample_cols <- grep("^(Zim|Zam)", names(tpm), value = TRUE)
if (length(sample_cols) < 4) stop("Fewer than four Zim/Zam samples were found.")

abundance <- tpm %>%
  select(all_of(sample_cols)) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ coalesce(.x, 0))) %>%
  as.matrix() %>%
  t()
colnames(abundance) <- tpm$Contig
rownames(abundance) <- sample_cols
relative <- sweep(abundance, 1, rowSums(abundance), "/")
group <- factor(ifelse(grepl("^Zim", rownames(relative)), "Zim-R", "Zam-R"))

alpha <- data.frame(
  Sample = rownames(relative),
  Group = group,
  Shannon = diversity(relative, index = "shannon")
)
write.csv(alpha, file.path(output_dir, "shannon_diversity.csv"), row.names = FALSE)

test <- wilcox.test(Shannon ~ Group, data = alpha)
p_alpha <- ggplot(alpha, aes(Group, Shannon, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, colour = "black") +
  geom_jitter(width = 0.12, size = 2) +
  annotate("text", x = 1.5, y = max(alpha$Shannon) * 1.05,
           label = paste0("Wilcoxon p = ", signif(test$p.value, 3))) +
  scale_fill_manual(values = c("Zim-R" = "#4E79A7", "Zam-R" = "#59A14F")) +
  labs(x = NULL, y = "Shannon diversity") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(output_dir, "shannon_diversity.png"), p_alpha,
       width = 5, height = 4, dpi = 300)

bray <- vegdist(relative, method = "bray")
permanova <- adonis2(bray ~ group, permutations = 999)
write.csv(as.data.frame(permanova), file.path(output_dir, "bray_curtis_permanova.csv"))

set.seed(42)
embedding <- Rtsne(as.matrix(bray), is_distance = TRUE,
                   dims = 2, perplexity = max(2, min(5, floor((nrow(relative) - 1) / 3))),
                   theta = 0)$Y
beta_plot <- data.frame(
  Sample = rownames(relative), Group = group,
  Axis1 = embedding[, 1], Axis2 = embedding[, 2]
)
p_beta <- ggplot(beta_plot, aes(Axis1, Axis2, colour = Group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Zim-R" = "#4E79A7", "Zam-R" = "#59A14F")) +
  labs(x = "t-SNE 1", y = "t-SNE 2", colour = NULL) +
  theme_classic()
ggsave(file.path(output_dir, "bray_curtis_tsne.png"), p_beta,
       width = 5, height = 4, dpi = 300)
