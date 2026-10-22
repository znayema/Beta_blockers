

#######################################################################
##APO miRNA Differential Expression Pipeline
## edgeR + Visualization (Boxplot, Paired Line, Heatmap)
## Author: Nayema
## Date: 04/28/2026
#######################################################################

rm(list = ls())

#######################################################################
# SECTION 1: LIBRARIES
#######################################################################

library(edgeR)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(haven)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(magick)
library(knitr)

#######################################################################
# SECTION 2: USER-DEFINED miRNA LIST
#######################################################################

mir_list <- c(
  "miR-100-5p",   "miR-106b-5p",  "miR-1285-5p",  "miR-132-3p",
  "miR-136-3p",   "miR-1538",     "miR-15b-5p",   "miR-194-3p",
  "miR-19a-3p",   "miR-20a-3p",   "miR-210-3p",   "miR-3199",
  "miR-338-3p",   "miR-374a-5p",  "miR-378c",     "miR-421",
  "miR-425-5p",   "miR-429",      "miR-4454",     "miR-452-5p",
  "miR-4654",     "miR-4667-5p",  "miR-4753-5p",  "miR-4999-5p",
  "miR-542-3p",   "miR-548ay-5p", "miR-548d-5p",  "miR-548n",
  "miR-580-3p",   "miR-590-3p",   "miR-6802-5p",  "miR-6842-5p",
  "miR-6859-5p",  "miR-6862-5p",  "miR-6882-5p",  "miR-6894-5p",
  "miR-873-5p",   "miR-874-5p",   "miR-885-5p"
)

#######################################################################
# SECTION 3: FILE PATHS (edit here only)
#######################################################################

PATH_EXPR    <- "/projects/larylab/data/beta_blocker_data/mirna/merged_mature_miRNA_expression.xlsx"
PATH_SAMPLES <- "/projects/larylab/data/beta_blocker_data/mirna/Sample List.xlsx"
PATH_TX      <- "/projects/larylab/data/beta_blocker_data/treatment_assignments.sas7bdat"

#######################################################################
# SECTION 4: LOAD RAW DATA
#######################################################################

expr_raw    <- read_excel(PATH_EXPR)
sample_raw  <- read_excel(PATH_SAMPLES)
txgroup_raw <- read_sas(PATH_TX)

#######################################################################
# SECTION 5: CLEAN METADATA
#######################################################################

# Clean treatment group — label as Atenolol / Placebo
txgroup <- txgroup_raw %>%
  mutate(
    ID    = as.character(ID),
    Group = ifelse(TXGROUP == 1, "Atenolol", "Placebo")
  )

# Clean sample list — standardize IDs and timepoint labels
sample_all <- sample_raw %>%
  mutate(
    RNAseq_ID = gsub("_R", "", `RNAseq ID`),
    Subject   = as.character(Subject),
    Timepoint = case_when(
      Timepoint %in% c("0", "Baseline", "baseline", "base")    ~ "Baseline",
      Timepoint %in% c("12", "Endpoint", "12 months")  ~ "Endpoint",
      TRUE                                               ~ NA_character_
    )
  ) %>%
  filter(!is.na(Timepoint)) %>%
  left_join(txgroup %>% select(ID, Group), by = c("Subject" = "ID")) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("Baseline", "Endpoint")),
    Group     = factor(Group,     levels = c("Atenolol", "Placebo"))
  )

cat("Samples after cleaning:\n")
print(table(sample_all$Group, sample_all$Timepoint))

#######################################################################
# SECTION 6: BUILD EXPRESSION MATRIX (all timepoints, mir_list only)
#######################################################################

# Build matrix, strip hsa- prefix, keep only mir_list miRNAs
expr_full <- expr_raw %>%
  mutate(miRNA_clean = sub("^hsa-", "", Mature.miRNA)) %>%
  filter(miRNA_clean %in% mir_list) %>%
  select(-Mature.miRNA) %>%
  column_to_rownames("miRNA_clean") %>%
  as.matrix()

cat("\nmiRNAs found from mir_list:", nrow(expr_full), "of", length(mir_list), "\n")

# Match to samples present in metadata
common_ids <- intersect(colnames(expr_full), sample_all$RNAseq_ID)
expr_full  <- expr_full[, common_ids, drop = FALSE]

sample_all <- sample_all %>%
  filter(RNAseq_ID %in% common_ids) %>%
  arrange(match(RNAseq_ID, common_ids))

cat("Samples matched to expression data:", length(common_ids), "\n")

#######################################################################
# SECTION 7: DIFFERENTIAL EXPRESSION — 12-MONTH ENDPOINT ONLY (edgeR)
#######################################################################

# Subset to endpoint samples only
sample_12mo <- sample_all %>% filter(Timepoint == "Endpoint")
expr_12mo   <- expr_full[, sample_12mo$RNAseq_ID, drop = FALSE]

# Build edgeR object
dge <- DGEList(
  counts = expr_12mo,
  group  = sample_12mo$Group
)

# Filter low-expression miRNAs
keep <- rowSums(cpm(dge) > 1) >= 3
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("\nmiRNAs retained after expression filter:", nrow(dge), "\n")

# Normalize and estimate dispersion
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# Exact test: Atenolol (1) vs Placebo (0)
de <- exactTest(dge, pair = c("Placebo", "Atenolol"))

# Compile results with FDR
results <- topTags(de, n = nrow(dge))$table %>%
  as.data.frame() %>%
  rownames_to_column("miRNA") %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>%
  arrange(PValue)

sig_results <- results %>% filter(FDR < 0.05)

cat("\nSignificant miRNAs (FDR < 0.05):", nrow(sig_results), "\n")
print(results)
print(sig_results)

# Save DE results
write.csv(results,     "/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/APO_DE_results_all.csv",         row.names = FALSE)
write.csv(sig_results, "/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/APO_DE_results_significant.csv",  row.names = FALSE)
cat("DE results saved.\n")

#######################################################################
# SECTION 8: NORMALIZE FULL DATA FOR PLOTTING (logCPM, both timepoints)
#######################################################################

dge_full    <- DGEList(counts = expr_full)
dge_full    <- calcNormFactors(dge_full)
logCPM_full <- cpm(dge_full, log = TRUE)

# Build long data frame for plotting
plot_df <- as.data.frame(logCPM_full) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "RNAseq_ID", values_to = "logCPM") %>%
  left_join(
    sample_all %>% select(RNAseq_ID, Subject, Timepoint, Group),
    by = "RNAseq_ID"
  ) %>%
  filter(!is.na(Group))

# Build long data frame with raw counts (for table + paired plot)
count_df <- as.data.frame(expr_full) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "RNAseq_ID", values_to = "counts") %>%
  left_join(
    sample_all %>% select(RNAseq_ID, Subject, Timepoint, Group),
    by = "RNAseq_ID"
  ) %>%
  filter(!is.na(Group))

#######################################################################
# SECTION 9: SUMMARY TABLE — Mean ± SD by Group × Timepoint (raw counts)
#######################################################################

#######################################################################
# SECTION 9: SUMMARY TABLE — Mean ± SD by Group × Timepoint (raw counts)
#######################################################################

# Verify both timepoints exist before building table
cat("Timepoints in count_df:\n")
print(table(count_df$Group, count_df$Timepoint, useNA = "ifany"))

summary_table <- count_df %>%
  # Ensure factors are characters for clean column naming
  mutate(
    Group     = as.character(Group),
    Timepoint = as.character(Timepoint)
  ) %>%
  group_by(miRNA, Group, Timepoint) %>%
  summarise(
    N    = n(),
    Mean = round(mean(counts, na.rm = TRUE), 2),
    SD   = round(sd(counts,   na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_SD  = paste0(Mean, " ± ", SD),
    col_name = paste(Group, Timepoint, sep = " | ")
  ) %>%
  select(miRNA, col_name, Mean_SD) %>%
  pivot_wider(
    names_from  = col_name,
    values_from = Mean_SD
  ) %>%
  # Reorder columns explicitly
  select(
    miRNA,
    any_of(c(
      "Atenolol | Baseline",
      "Atenolol | Endpoint",
      "Placebo | Baseline",
      "Placebo | Endpoint"
    ))
  )

print(kable(summary_table, format = "simple",
            caption = "miRNA Raw Counts: Mean ± SD by Group and Timepoint"))

write.csv(summary_table, "/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/APO_mirna_summary_mean_sd.csv", row.names = FALSE)
cat("Summary table saved: mirna_summary_mean_sd.csv\n")



table(sample_all$Group, sample_all$Timepoint, useNA = "ifany")
#######################################################################
# SECTION 10: COLOR PALETTES
#######################################################################mi

group_colors <- c(
  "Atenolol" = "#E07070",   # coral red
  "Placebo"  = "#5BC8C8"    # teal
)

#######################################################################
# SECTION 11: BOXPLOT FUNCTION (logCPM)
#######################################################################

plot_mirna_box <- function(mirna_name, df, colors) {
  df %>%
    filter(miRNA == mirna_name) %>%
    ggplot(aes(x = Timepoint, y = logCPM, fill = Group, color = Group)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha         = 0.7,
      position      = position_dodge(width = 0.75),
      width         = 0.6
    ) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      size     = 1.8,
      alpha    = 0.9
    ) +
    scale_fill_manual(values  = colors) +
    scale_color_manual(values = colors) +
    labs(
      title = mirna_name,
      y     = "logCPM",
      x     = "Timepoint",
      fill  = "Group",
      color = "Group"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10),
      plot.tag         = element_text(size = 14, face = "bold"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

#######################################################################
# SECTION 12: PAIRED LINE PLOT FUNCTION (raw counts, one line/patient)
#######################################################################

plot_mirna_paired <- function(mirna_name, df, colors) {
  
  df_sub <- df %>%
    filter(miRNA == mirna_name) %>%
    group_by(Subject) %>%
    filter(n_distinct(as.character(Timepoint)) == 2) %>%   # paired only
    ungroup()
  
  cat("Paired subjects for", mirna_name, ":",
      n_distinct(df_sub$Subject), "\n")
  cat("  Atenolol:", n_distinct(df_sub$Subject[df_sub$Group == "Atenolol"]), "\n")
  cat("  Placebo: ", n_distinct(df_sub$Subject[df_sub$Group == "Placebo"]),  "\n")
  
  ggplot(df_sub,
         aes(x = Timepoint, y = counts,
             group = Subject, color = Group)) +
    geom_line(alpha = 0.65, linewidth = 0.65) +
    geom_point(size = 2.2, alpha = 0.95) +
    scale_color_manual(values = colors) +
    scale_x_discrete(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = mirna_name,
      x     = "Timepoint",
      y     = "Counts",
      color = "Treatment"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position  = "right",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11),
      panel.background = element_rect(fill = "grey92", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

#######################################################################
# SECTION 13: GENERATE & SAVE BOXPLOTS (Panel A + B)
#######################################################################

mirna_A <- "miR-4667-5p"   # <-- change as needed
mirna_B <- "miR-19a-3p"    # <-- change as needed

pA <- plot_mirna_box(mirna_A, plot_df, group_colors) + labs(tag = "A")
pB <- plot_mirna_box(mirna_B, plot_df, group_colors) + labs(tag = "B")

pAB <- plot_grid(pA, pB, ncol = 1, align = "v")

ggsave("/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/figure4_boxplots.png", plot = pB,
       width = 4, height = 4, dpi = 300, bg = "white")
cat("Boxplots saved: figure4_boxplots.png\n")

#######################################################################
# SECTION 14: GENERATE & SAVE PAIRED LINE PLOT
#######################################################################

p_paired <- plot_mirna_paired("miR-19a-3p", count_df, group_colors)

ggsave("/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/figure_miR19a3p_paired.png", plot = p_paired,
       width = 5, height = 5, dpi = 300, bg = "white")
cat("Paired line plot saved: figure_miR19a3p_paired.png\n")

#######################################################################
# SECTION 15: HEATMAP — Mean logCPM, Z-score scaled rows
#######################################################################

# Build mean logCPM matrix: rows = miRNA, cols = Group × Timepoint
heatmap_mat <- plot_df %>%
  group_by(miRNA, Group, Timepoint) %>%
  summarise(mean_logCPM = mean(logCPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(col_label = paste(as.character(Group),
                           as.character(Timepoint), sep = " ")) %>%
  select(miRNA, col_label, mean_logCPM) %>%
  pivot_wider(names_from = col_label, values_from = mean_logCPM) %>%
  column_to_rownames("miRNA") %>%
  as.matrix()

# Order columns
col_order   <- c("Atenolol Baseline", "Atenolol Endpoint",
                 "Placebo Baseline",  "Placebo Endpoint")
col_order   <- col_order[col_order %in% colnames(heatmap_mat)]
heatmap_mat <- heatmap_mat[, col_order, drop = FALSE]

# Remove rows with NA or zero variance
heatmap_mat <- heatmap_mat[complete.cases(heatmap_mat), , drop = FALSE]
row_vars    <- apply(heatmap_mat, 1, var, na.rm = TRUE)
heatmap_mat <- heatmap_mat[!is.na(row_vars) & row_vars > 0, , drop = FALSE]
cat("miRNAs in heatmap:", nrow(heatmap_mat), "\n")

# Z-score scale rows
heatmap_scaled <- t(scale(t(heatmap_mat)))
cat("NA/Inf after scaling:",
    sum(is.na(heatmap_scaled) | is.infinite(heatmap_scaled)), "\n")

# Column annotations
annotation_col <- data.frame(
  Group     = ifelse(grepl("Atenolol", col_order), "Atenolol", "Placebo"),
  Timepoint = ifelse(grepl("Baseline", col_order), "Baseline", "Endpoint"),
  row.names = col_order
)

annotation_colors <- list(
  Group     = c(Atenolol = "#E07070", Placebo = "#5BC8C8"),
  Timepoint = c(Baseline = "#D9D9D9", Endpoint = "#555555")
)

pheatmap(
  heatmap_scaled,
  color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row      = 8,
  fontsize_col      = 9,
  angle_col         = 45,
  border_color      = NA,
  main              = "",
  filename          = "/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/figure4_heatmap.png",
  width             = 6,
  height            = 10,
  res               = 300
)
cat("Heatmap saved: figure4_heatmap.png\n")

#######################################################################
# SECTION 16: COMBINE BOXPLOTS + HEATMAP INTO ONE FIGURE
#######################################################################

img_ab <- image_read("figure4_boxplots.png")
img_c  <- image_read("figure4_heatmap.png")

target_h <- max(image_info(img_ab)$height, image_info(img_c)$height)

image_write(
  image_append(c(
    image_resize(img_ab, paste0("x", target_h)),
    image_resize(img_c,  paste0("x", target_h))
  )),
  "figure4_combined.png",
  quality = 95
)

data<-read_csv("/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/APO_mirna_summary_mean_sd.csv")
head(data)


library(tidyverse)

results <- count_df %>%
  group_by(miRNA, Group) %>%
  # Pivot to wide so Baseline and Endpoint are columns per subject
  pivot_wider(
    id_cols = c(miRNA, Subject, Group),
    names_from = Timepoint,
    values_from = counts
  ) %>%
  # Remove subjects missing either timepoint
  filter(!is.na(Baseline), !is.na(Endpoint)) %>%
  group_by(miRNA, Group) %>%
  summarise(
    n          = n(),
    mean_base  = mean(Baseline),
    sd_base    = sd(Baseline),
    mean_end   = mean(Endpoint),
    sd_end     = sd(Endpoint),
    p_value    = t.test(Baseline, Endpoint, paired = TRUE)$p.value,
    t_stat     = t.test(Baseline, Endpoint, paired = TRUE)$statistic,
    .groups    = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    across(c(mean_base, sd_base, mean_end, sd_end), ~round(.x, 2)),
    across(c(t_stat), ~round(.x, 3)),
    across(c(p_value, p_adj), ~signif(.x, 3)),
    `Baseline (Mean ± SD)`  = paste0(mean_base, " ± ", sd_base),
    `Endpoint (Mean ± SD)`  = paste0(mean_end,  " ± ", sd_end)
  ) %>%
  select(
    miRNA, Group, n,
    `Baseline (Mean ± SD)`,
    `Endpoint (Mean ± SD)`,
    t_stat, p_value, p_adj
  ) %>%
  arrange(Group, p_adj)

# View result
print(results)

# Optional: export to CSV
write_csv(results, "/projects/larylab/NAYEMA/Beta_blocker_project/mirna_beta_blockers/mirna/APO_mirna_paired_ttest_results.csv")
