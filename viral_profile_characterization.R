##################### Zim-Zam 
##### TPM normalization 
# Load necessary libraries
library(readxl)
library(tidyverse)

# 1. Load data
abund_raw <- read_excel("filtered_nonzero_contigs.xlsx")

# 2. Extract contig lengths from Contig names
abund_raw <- abund_raw %>%
  mutate(Length_bp = str_extract(Contig, "(?<=length_)[0-9]+") %>% as.numeric())

# 3. Keep only numeric abundance columns (assuming starts with "Zim"/"Zam")
abund_only <- abund_raw %>%
  select(Contig, Length_bp, starts_with("Zim"), starts_with("Zam"))

# 4. Compute TPM
compute_tpm <- function(df) {
  length_kb <- df$Length_bp / 1000  # Convert to kilobases
  counts <- df[, -c(1, 2)]          # Remove Contig and Length columns
  
  # Divide counts by length in kb to get RPK
  rpk <- sweep(counts, 1, length_kb, FUN = "/")
  
  # Compute scaling factor: sum of RPKs per sample / 1 million
  scaling_factors <- colSums(rpk, na.rm = TRUE) / 1e6
  
  # TPM = RPK / scaling factor
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")
  
  # Return TPM with Contig names
  tpm <- cbind(Contig = df$Contig, tpm)
  return(tpm)
}

# 5. Run TPM normalization
abund_tpm <- compute_tpm(abund_only)

# 6. Preview result
head(abund_tpm)

# 7. Save TPM data to CSV
write.csv(abund_tpm, file = "TPM_normalized_abundance.csv", row.names = FALSE)

######################### Diversity
#–– 0. Libraries ––#
library(readr)
library(dplyr)
library(writexl)

# 1. Read the TPM CSV file
df <- read_csv("TPM_normalized_abundance.csv")

# 2. Identify all abundance columns (those starting with "Zim" or "Zam")
ab_cols <- grep("^(Zim|Zam)", names(df), value = TRUE)

# 3. Filter out contigs that are zero in *all* of those columns
filtered <- df %>%
  filter(rowSums(across(all_of(ab_cols), ~ . != 0)) > 0)

# 4. Write the result to a new Excel file
write_xlsx(filtered, "TPM_filtered_nonzero_contigs.xlsx")

# 5. Optional: Print summary
cat(
  "✅ Kept", nrow(filtered), "of", nrow(df),
  "contigs (", nrow(df) - nrow(filtered),
  "all-zero contigs removed)\n"
)

######### Alpha diversity 
#–– Libraries ––#
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(patchwork)
library(grid)

#–– 1. Read TPM-normalized data for Shannon ––#
tpm_long <- read_csv("TPM_normalized_abundance.csv") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "abund") %>%
  mutate(group = ifelse(grepl("^Zim", sample), "Zim-R", "Zam-R")) %>%
  group_by(sample, group, Contig) %>%
  summarize(abund = sum(abund), .groups = "drop")

#–– 2. Read absolute counts data for Chao1 ––#
raw_counts <- read_excel("filtered_nonzero_contigs.xlsx") %>%
  select(Contig, starts_with("Zim"), starts_with("Zam"))

#–– 3. Build wide matrix for Shannon (TPM) ––#
tpm_wide <- tpm_long %>%
  pivot_wider(names_from = Contig, values_from = abund, values_fill = 0)

meta <- tpm_wide %>% select(sample, group)
mat_tpm <- tpm_wide %>% select(-sample, -group) %>% as.matrix()

#–– 4. Build wide matrix for Chao1 (raw counts) ––#
mat_counts <- raw_counts %>%
  select(-Contig) %>%
  t() %>%
  as.data.frame()

# Ensure sample names match
rownames(mat_counts) <- colnames(raw_counts)[-1]
mat_counts <- mat_counts[rownames(mat_counts) %in% meta$sample, ]

#–– 5. Compute alpha diversity ––#
alpha <- meta %>%
  mutate(
    Shannon = diversity(mat_tpm, "shannon"),
    Chao1 = estimateR(as.matrix(mat_counts))["S.chao1", meta$sample]
  )

#–– 6. Stats & labels ––#
p_shannon <- wilcox.test(Shannon ~ group, data = alpha)$p.value
p_chao1   <- wilcox.test(Chao1 ~ group, data = alpha)$p.value
lbl_shannon <- paste0("p = ", signif(p_shannon, 3))
lbl_chao1   <- paste0("p = ", signif(p_chao1, 3))

#–– 7. Sample labels ––#
alpha$group <- factor(alpha$group, levels = c("Zim-R", "Zam-R"))
ns <- alpha %>% count(group) %>% deframe()
xlabs <- c(
  paste0("Zim-R\n(", ns["Zim-R"], "/", ns["Zim-R"], ")"),
  paste0("Zam-R\n(", ns["Zam-R"], "/", ns["Zam-R"], ")")
)

#–– 8. Theme ––#
clean_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "black", fill = NA, size = 1),
    axis.title      = element_blank(),
    axis.text.x     = element_text(margin = margin(t = 6)),
    legend.position = "none"
  )

#–– 9. Plot function ––#
make_plot <- function(metric, p_label) {
  y_max <- max(alpha[[metric]], na.rm = TRUE) * 1.05
  ggplot(alpha, aes(x = group, y = .data[[metric]], fill = group)) +
    geom_boxplot(width = 0.6, color = "black", alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("Zim-R" = "#4E79A7", "Zam-R" = "#59A14F")) +
    scale_x_discrete(labels = xlabs) +
    annotate("text", x = 1.5, y = y_max, label = p_label, size = 4) +
    clean_theme +
    theme(axis.text = element_text(size = 12))
}

#–– 10. Generate plots ––#
p1 <- make_plot("Shannon", lbl_shannon)
p2 <- make_plot("Chao1", lbl_chao1)

#–– 11. Titles and layout ––#
title1 <- wrap_elements(grid::textGrob("Shannon Diversity", gp = gpar(fontface = "bold", fontsize = 16)))
title2 <- wrap_elements(grid::textGrob("Chao1 Richness", gp = gpar(fontface = "bold", fontsize = 16)))
titles <- title1 + title2 + plot_layout(ncol = 2)
final <- titles / (p1 + p2) + plot_layout(heights = c(0.1, 1))

#–– 12. Save and print ––#
ggsave("alpha_diversity_ZimZamR_final.png", final, width = 10, height = 6, dpi = 300)
print(final)


######### Beta diversity via t-SNE with stats annotation
#–– Libraries ––#
library(readr)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(Rtsne)
library(grid)

#–– 1. Load TPM-normalized abundance ––#
df_long <- read_csv("TPM_normalized_abundance.csv") %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "abund") %>%
  mutate(group = ifelse(grepl("^Zim", sample), "Zim-R", "Zam-R")) %>%
  group_by(sample, group, Contig) %>%
  summarize(abund = sum(abund), .groups = "drop")

#–– 2. Wide matrix + metadata ––#
df_wide <- df_long %>%
  pivot_wider(names_from = Contig, values_from = abund, values_fill = 0)

meta <- df_wide %>% select(sample, group)
mat  <- df_wide %>% select(-sample, -group) %>% as.matrix()
rownames(mat) <- df_wide$sample

#–– 3. Bray–Curtis distance (relative abundances) ––#
mat_rel   <- mat / rowSums(mat)
bray_dist <- vegdist(mat_rel, method = "bray")

#–– 4. PERMANOVA ––#
perm <- adonis2(bray_dist ~ group, data = meta, permutations = 999)
r2    <- perm$R2[1]
p_val <- perm$`Pr(>F)`[1]
stat_label <- paste0("R² = ", round(r2, 3), "\np = ", signif(p_val, 3))

#–– 5. t-SNE on Bray–Curtis distance ––#
set.seed(42)
tsne_out <- Rtsne(
  as.matrix(bray_dist),
  is_distance = TRUE,
  dims        = 2,
  perplexity  = 5,
  theta       = 0.0
)

#–– 6. Build t-SNE DataFrame ––#
tsne_df <- data.frame(
  sample = rownames(mat_rel),
  Dim1   = tsne_out$Y[, 1],
  Dim2   = tsne_out$Y[, 2]
) %>% left_join(meta, by = "sample")

#–– 7. Plot t-SNE with stats ––#
palette <- c("Zim-R" = "#4E79A7", "Zam-R" = "#59A14F")

p_tsne <- ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = palette) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = stat_label,
    hjust = -0.1, vjust = 1.1,
    size = 4
  ) +
  labs(
    title = "Beta Diversity (Bray–Curtis)",
    x = "tSNE1", y = "tSNE2",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text       = element_text(color = "black"),
    axis.ticks      = element_line(color = "black"),
    legend.position = "right"
  )

#–– 8. Save & Show ––#
ggsave("beta_tsne_ZimZamR_final.png", p_tsne, width = 6, height = 5, dpi = 300)
print(p_tsne)


######### Differential abundance via ANCOM-BC + Volcano plot
#!/usr/bin/env Rscript

# === USER CONFIGURATION =============================
input_xlsx <- "Zim-Zam.xlsx"
sheet_name <- "df_report"
output_dir <- "ancombc_results"
dir.create(output_dir, showWarnings = FALSE)

# === LOAD REQUIRED PACKAGES =========================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
  library(readr)
})

# === 1. READ CONTIG ABUNDANCE DATA ==================
df_report <- read_excel(input_xlsx, sheet = sheet_name)

# === 2. BUILD OTU MATRIX ============================
otu_df  <- df_report %>%
  select(Contig, starts_with("Zim"), starts_with("Zam")) %>%
  column_to_rownames("Contig")
otu_mat <- as.matrix(otu_df)
otu_tbl <- otu_table(otu_mat, taxa_are_rows = TRUE)

# === 3. CREATE SAMPLE METADATA ======================
samples <- colnames(otu_mat)
meta_df <- tibble(
  sample = samples,
  group  = ifelse(grepl("^Zim", samples), "Zim-R", "Zam-R")
) %>% column_to_rownames("sample")
sdata <- sample_data(meta_df)

# === 4. CREATE PHYLOSEQ OBJECT ======================
physeq <- phyloseq(otu_tbl, sdata)

# === 5. FILTER TAXA WITH <2 NON-ZERO COUNTS =========
physeq_filtered <- prune_taxa(rowSums(otu_table(physeq) > 0) >= 2, physeq)
cat("✔ Filtered to", ntaxa(physeq_filtered), "contigs with ≥2 non-zero values\n")

# === 6. RUN ANCOM-BC2 ===============================
res_ancom <- ancombc2(
  data         = physeq_filtered,
  fix_formula  = "group",
  p_adj_method = "BH",
  group        = "group",
  lib_cut      = 0,
  struc_zero   = FALSE,
  neg_lb       = FALSE,
  alpha        = 0.1
)

# === 7. TIDY RESULTS ================================
res_raw  <- as.data.frame(res_ancom$res) %>% rownames_to_column("Contig")
all_cols <- names(res_raw)

lfc_col <- setdiff(grep("^(beta_|lfc_|diff_)", all_cols, value = TRUE), grep("Intercept", all_cols, value = TRUE))[1]
qv_col  <- setdiff(grep("^(q_|padj|p_adj)", all_cols, value = TRUE), grep("Intercept", all_cols, value = TRUE))[1]

if (is.na(lfc_col) || is.na(qv_col)) {
  stop("❌ Could not auto-detect LFC or q-value columns.\nAvailable columns:\n", paste(all_cols, collapse = ", "))
} else {
  cat("✔ Using LFC column:   ", lfc_col, "\n")
  cat("✔ Using q-value col.: ", qv_col, "\n")
}

res_df <- res_raw %>%
  mutate(
    lfc        = .data[[lfc_col]],
    qval       = .data[[qv_col]],
    neglog10_q = -log10(qval),
    Significance = case_when(
      qval <  0.05 & lfc >  0 ~ "Enriched in Zam-R",
      qval <  0.05 & lfc <  0 ~ "Enriched in Zim-R",
      TRUE                    ~ "Not significant"
    )
  )

# === 8. SAVE RESULTS ================================
output_tsv <- file.path(output_dir, "ANCOMBC2_full_results.tsv")
write_tsv(res_df, output_tsv)
cat("✔ Results written to:", output_tsv, "\n")

# === 9. VOLCANO PLOT WITH LABELED TOP 5 PER GROUP ===========================
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)

# Identify top 5 most significant contigs per group
top_zim <- res_df %>%
  filter(Significance == "Enriched in Zim-R") %>%
  arrange(qval) %>%
  slice_head(n = 5)

top_zam <- res_df %>%
  filter(Significance == "Enriched in Zam-R") %>%
  arrange(qval) %>%
  slice_head(n = 5)

top_annotate <- bind_rows(top_zim, top_zam)

# Load metadata for VC, HostGenus, Lifestyle
meta_df <- read_excel("filtered_nonzero_contigs.xlsx")

# Format contig IDs to match metadata
top_annotate <- top_annotate %>%
  mutate(Contig = paste0("VC", Contig))

# Clean formatting to match metadata
top_annotate$Contig <- trimws(toupper(top_annotate$Contig))
meta_df <- meta_df %>%
  mutate(Contig = trimws(toupper(VC)))

# Generate annotation labels
top_labels <- meta_df %>%
  filter(Contig %in% top_annotate$Contig) %>%
  mutate(
    HostGenus = ifelse(is.na(vhost_hostGenus), "Unknown", vhost_hostGenus),
    RepMark   = ifelse(lifestyle == "temperate", "_T", "_V"),
    Label     = paste0(VC, "_", HostGenus, RepMark)
  ) %>%
  select(Contig, Label)

# Merge labels with top annotated data
top_annotate <- top_annotate %>%
  left_join(top_labels, by = "Contig")

# Define color palette
palette <- c(
  "Enriched in Zim-R" = "#4E79A7",
  "Enriched in Zam-R" = "#59A14F",
  "Not significant"   = "gray50"
)

# Generate volcano plot
volcano <- ggplot(res_df, aes(x = lfc, y = neglog10_q, color = Significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(alpha = 1, size = 3) +
  geom_text_repel(
    data = top_annotate,
    aes(label = Label),
    size = 2.6,
    fontface = "italic",
    color = "black",
    show.legend = FALSE,
    box.padding = 0.4,
    max.overlaps = 100,
    segment.size = 0.3,
    segment.color = "gray40"
  ) +
  scale_color_manual(values = palette, breaks = names(palette)) +
  labs(
    title = "ANCOM-BC2",
    x     = expression("Log"[2]*" Fold Change "),
    y     = expression("-Log"[10]*"(q-value)"),
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    legend.position  = "right"
  )

# Save the volcano plot
ggsave("ancombc_results/volcano_final_annotated_top10.png", plot = volcano, width = 7, height = 5.5, dpi = 300)

# Display in RStudio plot panel
print(volcano)


################################ Viral Contig–Metabolite Association Analysis

########## Step 0: Load Required Libraries ##########
library(readxl)
library(tidyverse)
library(pheatmap)
library(compositions)
library(reshape2)
library(ggplot2)
library(viridis)
library(writexl)

########## Step 1: Load Input Data ##########
abund_raw <- read_excel("filtered_nonzero_contigs.xlsx")
meta_raw  <- read_excel("Zim_Zam_SCFA_BA_v2.xlsx")

########## Step 2: Format Abundance Data ##########
abund <- abund_raw %>%
  select(Contig, starts_with("Zim"), starts_with("Zam")) %>%
  column_to_rownames("Contig")

########## Step 3: Format Metabolite Metadata ##########
meta <- meta_raw %>%
  column_to_rownames("id") %>%
  mutate(across(everything(), as.numeric))

########## Step 4: Align Shared Samples ##########
shared_samples <- intersect(colnames(abund), rownames(meta))
abund <- abund[, shared_samples]
meta  <- meta[shared_samples, ]

########## Step 5: Filter Contigs Present in ≥3 Samples ##########
abund <- abund[rowSums(abund > 0) >= 3, ]

########## Step 6: CLR Transformation on Viral Abundance ##########
abund_clr <- abund %>%
  as.matrix() %>%
  t() %>%
  apply(1, function(x) clr(x + 1e-5)) %>%
  t() %>%
  as.data.frame()
colnames(abund_clr) <- rownames(abund)
rownames(abund_clr) <- shared_samples

########## Step 7: Z-score Transformation on Metabolites ##########
meta_z <- meta %>%
  mutate(across(everything(), scale)) %>%
  as.data.frame()

########## Step 8: Spearman Correlation Matrix ##########
n_contigs <- ncol(abund_clr)
n_mets    <- ncol(meta_z)

cor_mat  <- matrix(NA, nrow = n_contigs, ncol = n_mets)
pval_mat <- matrix(NA, nrow = n_contigs, ncol = n_mets)
rownames(cor_mat) <- rownames(pval_mat) <- colnames(abund_clr)
colnames(cor_mat) <- colnames(pval_mat) <- colnames(meta_z)

for (i in seq_len(n_contigs)) {
  for (j in seq_len(n_mets)) {
    x <- abund_clr[, i]
    y <- meta_z[, j]
    if (sd(x) > 0 && sd(y) > 0) {
      test <- cor.test(x, y, method = "spearman")
      cor_mat[i, j]  <- test$estimate
      pval_mat[i, j] <- test$p.value
    }
  }
}

########## Step 9: Optional Heatmap of All Correlations ##########
filtered_cor <- cor_mat[rowSums(!is.na(cor_mat)) > 0, colSums(!is.na(cor_mat)) > 0]

if (nrow(filtered_cor) >= 2 && ncol(filtered_cor) >= 2) {
  pheatmap(filtered_cor,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           na_col = "grey90",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "All CLR-Transformed Viral Contigs vs Z-scored Metabolites\n(Spearman Correlation)")
}

########## Step 10: Save Matrices ##########
write.csv(cor_mat, "correlation_matrix_all.csv")
write.csv(pval_mat, "pvalues_matrix_all.csv")

########## Step 11: Process Correlation Data ##########
cor_mat[abs(cor_mat) < 0.4] <- NA
pval_mat[is.na(cor_mat)] <- NA

r_long    <- melt(as.matrix(cor_mat), varnames = c("Contig", "Metabolite"), value.name = "r")
pval_long <- melt(as.matrix(pval_mat), varnames = c("Contig", "Metabolite"), value.name = "p")

data_long <- left_join(r_long, pval_long, by = c("Contig", "Metabolite")) %>%
  mutate(
    p_adj = p.adjust(p, method = "fdr"),
    sig_star = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    r_label = ifelse(!is.na(r), sprintf("%.2f%s", r, sig_star), "")
  )

########## Step 12: Top 20 Contigs by Significant Associations ##########
top_contigs <- data_long %>%
  filter(!is.na(r)) %>%
  group_by(Contig) %>%
  summarise(n_strong = n()) %>%
  arrange(desc(n_strong)) %>%
  slice_head(n = 20) %>%
  pull(Contig)

data_top20 <- data_long %>%
  filter(Contig %in% top_contigs)

########## Step 13: Clean Labels for Top 20 ##########
meta_df <- read_excel("filtered_nonzero_contigs.xlsx")

contig_info <- meta_df %>%
  filter(Contig %in% top_contigs) %>%
  mutate(
    HostGenus = ifelse(is.na(vhost_hostGenus), "Unknown", vhost_hostGenus),
    RepMark   = ifelse(lifestyle == "temperate", "_T", ""),
    Label     = paste0(VC, "_", HostGenus, RepMark)
  ) %>%
  select(Contig, Label)

data_top20_labeled <- data_top20 %>%
  left_join(contig_info, by = "Contig") %>%
  mutate(
    Contig = Label,
    text_color = ifelse(r >= 0.4, "white", "black")
  )

####### Step 14: Final Plot – Top 20 Contigs Heatmap
charcoal_blue <- colorRampPalette(c("#f7f9fb", "#b0c4de", "#2f4f4f"))

ggplot(data_top20_labeled, aes(x = Metabolite, y = Contig, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = r_label, color = text_color), size = 2.5, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = charcoal_blue(100),
    na.value = "white",
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title  = element_text(hjust = 0.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.2, "cm"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8)
  ) +
  labs(
    title = "Top 20 Viral Contigs Correlation with Metabolites",
    x = "Metabolite",
    y = "Viral Contig"
  )

# === 9. VOLCANO PLOT WITH LABELED TOP 5 PER GROUP ===========================

# Identify top 5 most significant contigs per group
top_zim <- res_df %>%
  filter(Significance == "Enriched in Zim-R") %>%
  arrange(qval) %>%
  slice_head(n = 5)

top_zam <- res_df %>%
  filter(Significance == "Enriched in Zam-R") %>%
  arrange(qval) %>%
  slice_head(n = 5)

top_annotate <- bind_rows(top_zim, top_zam)

# Load metadata for VC, HostGenus, Lifestyle
meta_df <- read_excel("filtered_nonzero_contigs.xlsx")

top_labels <- meta_df %>%
  filter(Contig %in% top_annotate$Contig) %>%
  mutate(
    HostGenus = ifelse(is.na(vhost_hostGenus), "Unknown", vhost_hostGenus),
    RepMark   = ifelse(lifestyle == "temperate", "_T", "_V"),
    Label     = paste0(VC, "_", HostGenus, RepMark)
  ) %>%
  select(Contig, Label)

# Merge annotated labels with top 10 points
top_annotate <- top_annotate %>%
  left_join(top_labels, by = "Contig")

# Plotting
palette <- c(
  "Enriched in Zim-R" = "#4E79A7",
  "Enriched in Zam-R" = "#59A14F",
  "Not significant"   = "gray50"
)

volcano <- ggplot(res_df, aes(x = lfc, y = neglog10_q, color = Significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(alpha = 0.8, size = 2) +
  geom_text(
    data = top_annotate,
    aes(label = Label),
    size = 3,
    vjust = -1,
    fontface = "italic",
    color = "black",
    show.legend = FALSE
  ) +
  scale_color_manual(values = palette, breaks = names(palette)) +
  labs(
    title = "ANCOM-BC2",
    x     = expression("Log"[2]*" Fold Change "),
    y     = expression("-Log"[10]*"(q-value)"),
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    legend.position  = "right"
  )

############################# Replication cycle and host phyla
#–––– Libraries ––#
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#–––– Step 1: Load TPM and Metadata ––#
tpm  <- read_excel("TPM_filtered_nonzero_contigs.xlsx")
meta <- read_excel("filtered_nonzero_contigs.xlsx")

#–––– Step 2: Extract Zim and Zam sample names ––#
zim_samples <- grep("^Zim", colnames(tpm), value = TRUE)
zam_samples <- grep("^Zam", colnames(tpm), value = TRUE)
all_samples <- c(zim_samples, zam_samples)

#–––– Step 3: Merge TPM and metadata ––#
tpm_long <- tpm %>%
  select(Contig, all_of(all_samples)) %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "TPM")

merged <- tpm_long %>%
  left_join(meta, by = "Contig") %>%
  mutate(Group = ifelse(grepl("^Zim", Sample), "Zim_R", "Zam_R"))

#–––– Custom Palette ––#
palette <- c("Zim_R" = "#4E79A7", "Zam_R" = "#59A14F")

#–––– 1. Host Phylum Composition (Excluding Unknown) ––#
phylum_data <- merged %>%
  filter(!is.na(vhost_hostPhylum)) %>%  # Exclude Unknown
  group_by(Sample, Group, HostPhylum = vhost_hostPhylum) %>%
  summarise(RelAbund = sum(TPM), .groups = "drop")

ggplot(phylum_data, aes(x = HostPhylum, y = RelAbund, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.8) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.format", label.y = max(phylum_data$RelAbund) * 1.05) +
  scale_fill_manual(values = palette) +
  labs(
    title = "Host Phylum Composition",
    y = "Relative Abundance (%)",
    x = "Host Phylum"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  )

#–––– 2. Replication Cycle Composition ––#
lifestyle_data <- merged %>%
  filter(!is.na(lifestyle)) %>%
  group_by(Sample, Group, lifestyle) %>%
  summarise(RelAbund = sum(TPM), .groups = "drop")

ggplot(lifestyle_data, aes(x = lifestyle, y = RelAbund, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.6),
               outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.8) +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.format", label.y = max(lifestyle_data$RelAbund) * 1.05) +
  scale_fill_manual(values = palette) +
  labs(
    title = "Replication Cycle Composition",
    y = "Relative Abundance (%)",
    x = "Replication cycle"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  )


########################### Compare Contig Presence Using a Venn Diagram
#–––– Libraries ––––#
library(readxl)
library(dplyr)
library(VennDiagram)
library(grid)

#–––– Step 1: Load TPM Data ––––#
tpm <- read_excel("TPM_filtered_nonzero_contigs.xlsx")

#–––– Step 2: Extract Zim and Zam Samples ––––#
zim_samples <- grep("^Zim", colnames(tpm), value = TRUE)
zam_samples <- grep("^Zam", colnames(tpm), value = TRUE)

#–––– Step 3: Determine Presence per Group ––––#
contig_presence <- tpm %>%
  mutate(
    Zim_present = rowSums(across(all_of(zim_samples)) > 0) > 0,
    Zam_present = rowSums(across(all_of(zam_samples)) > 0) > 0
  ) %>%
  select(Contig, Zim_present, Zam_present)

#–––– Step 4: Create Sets ––––#
zim_set <- contig_presence$Contig[contig_presence$Zim_present]
zam_set <- contig_presence$Contig[contig_presence$Zam_present]
intersect_set <- intersect(zim_set, zam_set)

#–––– Step 5: Draw Venn Diagram ––––#
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = length(zim_set),
  area2 = length(zam_set),
  cross.area = length(intersect_set),
  category = c("", ""),                      # disable default labels
  fill = c("#4E79A7", "#59A14F"),            # Zim and Zam fill colors
  lty = "solid",                             # border line type
  lwd = 4,                                   # border line width
  col = c("white", "white"),                 # circle outline color
  alpha = c(0.6, 0.6),                       # transparency of fill
  cex = 2.2,                                 # size of numbers inside
  fontfamily = "sans",                       # font
  euler.d = FALSE,                           # perfect circles
  scaled = FALSE                             # disable size scaling
)

#–––– Step 6: Add manual group labels above each circle ––––#
# You can slightly tweak x/y values for positioning
grid.text("Zim-R",
          x = unit(0.32, "npc"),
          y = unit(0.9, "npc"),
          gp = gpar(col = "gray30", fontsize = 16))

grid.text("Zam-R",
          x = unit(0.68, "npc"),
          y = unit(0.9, "npc"),
          gp = gpar(col = "gray30", fontsize = 16))


####################################### Network analyses Spearman with bootstraping and FDR
# Load libraries
library(readr)
library(readxl)
library(dplyr)
library(igraph)
library(ggraph)
library(scales)
library(ggplot2)

# Load correlation results
df <- read.delim("/Users/ali/Desktop/Zim-Zam/spearman_two_layer_bootstrap_fdr_rho0.3.tsv")

# Filter based on correlation strength and reliability
df_filtered <- df %>%
  filter(SpearmanRho >= 0.5, FDR <= 0.05, BootstrapPassRatio >= 0.95)

# Identify node types (plural)
get_group <- function(x) {
  if (grepl("^b_", x)) return("Bacteria")
  if (grepl("^v_", x)) return("Viruses")
  if (grepl("^m_", x)) return("Metabolites")
  return("Other")
}

# Keep only cross-domain edges
edges_cross <- df_filtered %>%
  mutate(group1 = sapply(Feature1, get_group),
         group2 = sapply(Feature2, get_group)) %>%
  filter(group1 != group2) %>%
  select(from = Feature1, to = Feature2, weight = SpearmanRho)

# Build graph
g <- graph_from_data_frame(edges_cross, directed = FALSE)
V(g)$group <- sapply(V(g)$name, get_group)

# Extract largest connected component with all three groups
comps <- components(g)
valid_comps <- which(sapply(unique(comps$membership), function(cid) {
  nodes <- names(comps$membership[comps$membership == cid])
  groups <- unique(V(g)$group[V(g)$name %in% nodes])
  all(c("Bacteria", "Viruses", "Metabolites") %in% groups)
}))
tripartite_nodes <- names(comps$membership[comps$membership %in% valid_comps])
tripartite_subgraph <- induced_subgraph(g, tripartite_nodes)
largest <- which.max(sizes(components(tripartite_subgraph)))
g_tri <- induced_subgraph(tripartite_subgraph, V(tripartite_subgraph)[components(tripartite_subgraph)$membership == largest])

# Load mapping files
bacteria_map_raw <- read_csv("/Users/ali/Desktop/Zim-Zam/OTUs_renamed_step1.csv", show_col_types = FALSE)
virus_map <- read_xlsx("/Users/ali/Desktop/Zim-Zam/filtered_nonzero_contigs.xlsx")

# Bacteria label extraction
extract_clean_label <- function(tax_string, fallback_id) {
  tax <- trimws(strsplit(tax_string, ";")[[1]])
  names(tax) <- c("k", "p", "c", "o", "f", "g", "s")[seq_along(tax)]
  bad_keywords <- c("Unclassified", "uncultured", "unknown", "metagenome", "group", "UCG", "clade", "subgroup", "NK4A136", "R-7", "RC9 gut", "NK4A214")
  for (lvl in c("g", "f", "o", "c", "p")) {
    val <- tax[lvl]
    if (!is.na(val) && val != "") {
      val_clean <- gsub(".*__", "", val)
      val_clean <- gsub(" [^;]*", "", val_clean)
      if (val_clean != "" && !any(grepl(paste(bad_keywords, collapse = "|"), val_clean, ignore.case = TRUE))) {
        return(paste0(lvl, "_", val_clean))
      }
    }
  }
  return(fallback_id)
}

bacteria_labels <- bacteria_map_raw %>%
  rename(FeatureID = `##OTU ID`) %>%
  rowwise() %>%
  mutate(
    clean_label = extract_clean_label(taxonomy, FeatureID),
    b_id = paste0("b_", FeatureID)
  ) %>%
  ungroup() %>%
  select(b_id, clean_label)

# Virus label cleaning
virus_map <- virus_map %>%
  mutate(v_id = paste0("v_", Contig)) %>%
  rowwise() %>%
  mutate(
    host_taxa = coalesce(
      vhost_hostSpecies, vhost_hostGenus, vhost_hostFamily,
      vhost_hostOrder, vhost_hostPhylum, vhost_hostName, "UnknownHost"
    ),
    lifestyle_short = case_when(
      lifestyle == "temperate" ~ "T",
      lifestyle == "virulent" ~ "V",
      TRUE ~ "X"
    ),
    new_label = gsub("^VCVC", "VC", paste0("VC", VC, "_", host_taxa, "_", lifestyle_short))
  ) %>%
  ungroup() %>%
  select(v_id, new_label)

# Relabel nodes
deg <- degree(g_tri)
rename_df <- data.frame(original = V(g_tri)$name, stringsAsFactors = FALSE) %>%
  mutate(group = sapply(original, get_group)) %>%
  mutate(m_label = ifelse(group == "Metabolites", gsub("^m_", "", original), NA)) %>%
  left_join(bacteria_labels, by = c("original" = "b_id")) %>%
  left_join(virus_map, by = c("original" = "v_id")) %>%
  mutate(label = case_when(
    !is.na(new_label) ~ new_label,
    !is.na(clean_label) ~ clean_label,
    !is.na(m_label) ~ m_label,
    TRUE ~ NA_character_
  ))

V(g_tri)$label <- rename_df$label

# Degree-based sizing: 3 levels
V(g_tri)$size_group <- cut(deg,
                           breaks = c(-Inf, 4, 8, Inf),
                           labels = c("≤4", "5–8", ">8")
)
V(g_tri)$degree <- deg

# Final plot
set.seed(123)
ggraph(g_tri, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.5, color = "gray50") +
  geom_node_point(aes(color = group, size = size_group)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 2.5) +
  scale_size_manual(
    name = "Node Centrality",
    values = c("≤4" = 2, "5–8" = 4, ">8" = 6)
  ) +
  scale_color_manual(
    name = "Node Type",
    values = c(
      "Bacteria"    = "#688D74",
      "Metabolites" = "#E3B187",
      "Viruses"     = "#683A33"
    )
  ) +
  scale_edge_width(range = c(0.3, 1.2), name = "Spearman ρ") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


########################### ML of the effect of viruses on bacterial metabolites 
# Load libraries
library(ggplot2)
library(readr)

# Load data
df <- read_csv("virus_effects_on_metabolites.csv")

# Fix ordering for the plot
df$Metabolite <- factor(df$Metabolite, levels = df$Metabolite[order(df$Delta_R2_Virus_on_Bacteria)])

# Plot
ggplot(df, aes(x = Metabolite, y = Delta_R2_Virus_on_Bacteria, fill = EffectType)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -0.05, linetype = "dashed", color = "darkred", size = 0.5) +
  scale_fill_manual(values = c("Additive" = "darkgreen", "Neutral" = "gray", "Suppressive" = "red")) +
  labs(
    title = "Effect of Adding Viral Features to Bacterial Models",
    y = expression(Delta~R^2~"(Combined − Bacteria Only)"),
    x = "Metabolite"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))