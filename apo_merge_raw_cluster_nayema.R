#######################################################################
##APO miRNA Differential Expression Pipeline
## edgeR + Visualization (Boxplot, Paired Line, Heatmap)
## + BMD MANOVA extension (outlier removal + additional BMD traits)
## Author: Nayema
## Date: 04/28/2026
#######################################################################

rm(list = ls())
renv::activate()
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
library(readr)
library(tidyverse)
library(lme4)
library(factoextra)
library(car)

#######################################################################
# SECTION 3: FILE PATHS (edit here only)
#######################################################################

PATH_EXPR    <- "/projects/larylab/data/beta_blocker_data/mirna/merged_mature_miRNA_expression.xlsx"
PATH_SAMPLES <- "/projects/larylab/data/beta_blocker_data/mirna/Sample List.xlsx"
PATH_TX      <- "/projects/larylab/data/beta_blocker_data/treatment_assignments.sas7bdat"
PATH_CLIN    <- "/projects/larylab/data/beta_blocker_data/DXA_data/apodxa_corrected_20260309.sas7bdat"
PATH_BTM     <- "/projects/larylab/data/beta_blocker_data/1-7-26 APO CTx data.xlsx"
PATH_OSTEOCLAST <- "/projects/larylab/data/beta_blocker_data/Copy of maine_geno_pheno_updated.xlsx"
PATH_GENO       <- "/projects/larylab/data/beta_blocker_data/four_snps_geno.rds"

#######################################################################
# SECTION 3B: CONFIG — outlier samples + additional BMD traits
#######################################################################

EXCLUDE_RNASEQ_IDS <- c("DM132", "DM142", "DM17")


BMD_TRAIT_COLS <- c(
  HTOTBMD  = "hipbmd",     # total hip BMD
  STOTBMD  = "spinebmd",   # total spine BMD
  RTOTBMD  = "rtotbmd",    # radius total BMD
  UTOTBMD  = "utotbmd",    # ulna total BMD
  RUTOTBMD = "rutotbmd"    # radius+ulna total BMD
)

# Human-readable labels for figure/table titles
BMD_TRAIT_LABELS <- c(
  hipbmd    = "Hip Total BMD",
  spinebmd  = "Spine Total BMD",
  rtotbmd   = "Radius Total BMD",
  utotbmd   = "Ulna Total BMD",
  rutotbmd  = "Radius+Ulna Total BMD"
)

# Where figures/tables get written
OUT_DIR <- "/projects/larylab/beta_blockers/mirna/results/bmd_outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

#######################################################################
# SECTION 4: LOAD RAW DATA
#######################################################################

expr        <- read_excel(PATH_EXPR)
sample  <- read_excel(PATH_SAMPLES)
txgroup <- read_sas(PATH_TX)
clin <- read_sas(PATH_CLIN)
btm <- read_excel(PATH_BTM)
osteoclast <- read_excel(PATH_OSTEOCLAST)
geno <- readRDS(PATH_GENO)

#######################################################################
# SECTION 4B: EXCLUDE OUTLIER SAMPLES (map RNAseq ID -> Subject/APO id)
#######################################################################

exclude_subjects <- sample %>%
  filter(`RNAseq ID` %in% EXCLUDE_RNASEQ_IDS) %>%
  pull(Subject) %>%
  unique()

cat("Excluding subjects (mapped from", paste(EXCLUDE_RNASEQ_IDS, collapse = ", "),
    "):", paste(exclude_subjects, collapse = ", "), "\n")

# Sanity check: make sure every requested RNAseq ID actually resolved to a subject
missing_map <- setdiff(EXCLUDE_RNASEQ_IDS, sample$`RNAseq ID`)
if (length(missing_map) > 0) {
  warning("These RNAseq IDs were not found in `sample` and were not excluded: ",
          paste(missing_map, collapse = ", "))
}

# Drop from every raw table that carries an id, before any joins/renames
clin       <- clin       %>% filter(!ID %in% exclude_subjects)
txgroup    <- txgroup    %>% filter(!ID %in% exclude_subjects)
btm        <- btm        %>% filter(!`Subject ID` %in% exclude_subjects)
osteoclast <- osteoclast %>% filter(!Subject.ID...5 %in% exclude_subjects) # pre-"APO"-prefix; see rename below
geno       <- geno       %>% filter(!Subject.ID %in% exclude_subjects)

#######################################################################
# SECTION 5: CLEAN METADATA
#######################################################################

# subset CTX + additional BMD trait data
clin <- clin %>%
  select(`ID`, all_of(names(BMD_TRAIT_COLS)), `SITE`, `VISITN`) %>%
  rename(id = `ID`, site = `SITE`, time = `VISITN`) %>%
  rename(!!!setNames(names(BMD_TRAIT_COLS), unname(BMD_TRAIT_COLS)))

btm <- btm %>%
  select(`Subject ID`,`Site name`,`visit month`,`Sample Data Display Value`) %>%
  rename(id=`Subject ID`,site=`Site name`,time=`visit month`,ctx=`Sample Data Display Value`)

osteoclast <- read_excel(PATH_OSTEOCLAST)
osteoclast <- osteoclast %>%
  filter(!Subject.ID...5 %in% exclude_subjects) %>%
  rename(id = `Subject.ID...5`) %>%
  rename(CTX_veh = "Vehicle", CTX_aten = `Aten 15`) %>%
  rename(adrb2 = "adrb2_I_think") %>%
  select(id,CTX_veh,CTX_aten,rs1801252_a_g,rs12414657_t_c,adrb2)

osteoclast <- osteoclast %>%
  mutate(id = paste("APO",id,sep=""))

# Clean treatment group — label as Atenolol / Placebo
txgroup <- txgroup %>%
  mutate(
    ID    = as.character(ID),
    Group = ifelse(TXGROUP == 1, "Atenolol", "Placebo")
  )

btmwide <- btm %>%
  pivot_wider(id_cols=id,names_from = time,values_from = ctx) %>%
  rename(ctx0 = `0`,ctx6 = `6`, ctx12 = `12`, ctx18 = `18`, ctx24 = `24`) %>%
  mutate(ctx0 = as.numeric(ctx0), ctx6 = as.numeric(ctx6), ctx12 = as.numeric(ctx12),
         ctx18 = as.numeric(ctx18), ctx24 = as.numeric(ctx24)) %>%
  left_join(osteoclast,by="id") %>%
  left_join(txgroup,by=c("id" = "ID")) %>%
  select(-md2randid,-TXGROUP)

# plot timecourse of CTX by individual colored by group
btm_treat <- btmwide %>%
  dplyr::select(id,ctx0,ctx6,ctx12,ctx18,ctx24,Group) %>%
  pivot_longer(cols=c("ctx0","ctx6","ctx12","ctx18","ctx24"),names_to="visit",values_to="ctx") %>%
  mutate(ctx = as.numeric(ctx)) %>%
  mutate(t = as.numeric(gsub("ctx","",visit)))

# make geno wide
geno_wide <- geno %>%
  dplyr::select(Subject.ID,RSID,Effect_Allele_Count) %>%
  rename(id = Subject.ID) %>%
  rename(SNP = RSID) %>%
  rename(geno = Effect_Allele_Count) %>%
  pivot_wider(id_cols=id,names_from = SNP,values_from = geno)

#######################################################################
# SECTION: response/rebound classification (unchanged from original,
# just operating on data with the 3 outlier subjects already removed)
#######################################################################

first <- btm_treat %>% filter(t <= 6)
second <- btm_treat %>% filter(t > 6)
third <- btm_treat %>% filter(t %in% c(0,24))

ids = unique(first$id)
n = length(ids)
param = array(0,dim=c(n,3))
for (i in 1:n) {
  mod1 = lm(ctx~t,data=first[first$id==ids[i],])
  param[i,1] = mod1$coef[2]
  mod2 = lm(ctx~t,data=second[second$id==ids[i],])
  param[i,2] = mod2$coef[2]
  mod3 = lm(ctx~t,data=third[third$id==ids[i],])
  param[i,3] = mod3$coef[2]
}
paramdf = data.frame(param)
names(paramdf)=c("Response","Rebound","Total")
paramdf$id = ids
paramdf <- paramdf %>%
  left_join(txgroup,by=c("id" = "ID"))

paramdf <- paramdf %>%
  mutate(respond = (Response < -0.02) & (Rebound < 0.01)) %>%
  mutate(respond2 = (Response < -0.02))

btm_treat1 <- btm_treat %>%
  left_join(paramdf,by=c("id","Group")) %>%
  mutate(Group2 = ifelse(Group=="Placebo","Placebo",ifelse(respond,"Responder","NR"))) %>%
  mutate(Group3 = ifelse(Group=="Placebo","Placebo",ifelse(respond2,"Responder","NR"))) %>%
  right_join(geno_wide,by="id") %>%
  mutate(rs1042714_treat = ifelse(rs1042714==0,"GG",ifelse(rs1042714==1,"GC","CC"))) %>%
  mutate(Group4 = ifelse(Group=="Placebo","Placebo",rs1042714_treat))

# do manova to test if CTX trajectory is different by group
btm_treat1_wide <- btm_treat1 %>%
  mutate(rs1042714 = ifelse(rs1042714==0,"GG",ifelse(rs1042714==1,"GC","CC"))) %>%
  filter(!is.na(t)) %>%
  mutate(t = paste("ctx_",t,sep="")) %>%
  dplyr::select(id,Group,t,rs1042714,ctx) %>%
  pivot_wider(id_cols=c("id","Group","rs1042714"),names_from = "t",values_from = "ctx")

#######################################################################
# SECTION: MANOVA on CTX trajectory (original, unchanged in logic)
#######################################################################

ctx <- with(btm_treat1_wide, cbind(ctx_0,ctx_6,ctx_12,ctx_18,ctx_24))
fit <- manova(ctx ~ Group*rs1042714, data = btm_treat1_wide)
summary(fit, test = "Wilks")

mod <- lm(ctx ~ Group*rs1042714, data = btm_treat1_wide)
idata <- data.frame(Time = factor(1:5))
rm.manova <- Anova(mod, idata = idata, idesign = ~Time, type = "III")
summary(rm.manova, multivariate = TRUE)

#######################################################################
# SECTION: BMD ANALYSIS — mirrors the CTX analysis pipeline above,
# applied to all 5 "Total" BMD traits, on the sample with the 3 outlier
# subjects already excluded.
#
# ADAPTATION NOTE: the CTX pipeline defines "Response"/"Rebound" slopes
# (0-6mo / 6-24mo) to capture a biphasic drop-then-rebound in a fast-
# turnover blood marker, then classifies subjects as "responders" using
# CTX-specific cutoffs (-0.02 / 0.01 ng/mL per month). BMD is a slow-
# changing structural measure with no equivalent validated cutoff, so
# below we:
#   (a) DO compute the same three slopes (Response/Rebound/Total) per
#       subject per trait, for a like-for-like comparison structure,
#   (b) do NOT invent a BMD "responder" cutoff -- instead we test Group
#       differences directly on the continuous slopes (t-test + boxplot,
#       same tools she used for Response~Group on CTX),
#   (c) test genotype association via ANOVA on the continuous slopes
#       (aov(slope ~ genotype)) rather than a cutoff-based crosstab.
#######################################################################

# Base dataset: BMD traits joined to treatment group + rs1042714 genotype
clin_base <- clin %>%
  left_join(txgroup, by = c("id" = "ID")) %>%
  left_join(geno_wide, by = "id") %>%
  mutate(rs1042714_cat = ifelse(rs1042714 == 0, "GG",
                                ifelse(rs1042714 == 1, "GC", "CC")))

# Guard against duplicate id x time rows (e.g. bilateral hip/forearm scans)
# before pivoting wide -- average if duplicates exist.
clin_base <- clin_base %>%
  group_by(id, time, Group, rs1042714_cat) %>%
  summarize(across(all_of(unname(BMD_TRAIT_COLS)), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

# Long format across all BMD traits -- used for plotting, slopes, tables
clin_long <- clin_base %>%
  pivot_longer(cols = all_of(unname(BMD_TRAIT_COLS)),
               names_to = "trait", values_to = "bmd") %>%
  filter(!is.na(time), time %in% c(0, 6, 12, 18, 24)) %>%   # drop the stray "25" visit (n<15, likely an unscheduled extra visit)
  mutate(trait_label = BMD_TRAIT_LABELS[trait])

build_trait_wide <- function(data, trait, times = c(0,6,12,18,24)) {
  data %>%
    filter(time %in% times, !is.na(time)) %>%
    select(id, Group, rs1042714 = rs1042714_cat, time, value = all_of(trait)) %>%
    mutate(time = paste0(trait, "_", time)) %>%
    pivot_wider(id_cols = c(id, Group, rs1042714), names_from = time, values_from = value)
}

#######################################################################
#######################################################################

plot_trait_trajectory <- function(long_df, trait, label) {
  d <- long_df %>% filter(trait == !!trait)
  d_summary <- d %>%
    group_by(Group, time) %>%
    summarize(mean = mean(bmd, na.rm = TRUE),
              se   = sd(bmd, na.rm = TRUE) / sqrt(sum(!is.na(bmd))),
              .groups = "drop")
  
  p <- ggplot(d_summary, aes(x = time, y = mean, color = Group)) +
    geom_ribbon(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se, fill = Group),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(Atenolol = "#D55E00", Placebo = "#0072B2")) +
    scale_fill_manual(values  = c(Atenolol = "#D55E00", Placebo = "#0072B2")) +
    scale_x_continuous(breaks = c(0,6,12,18,24)) +
    labs(title = paste0(label, ": mean \u00b1 95% CI by treatment group"),
         x = "Visit month", y = expression(paste("Mean BMD (g/cm"^2, ")"))) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  
  ggsave(file.path(OUT_DIR, paste0("fig_trajectory_", trait, ".png")),
         plot = p, width = 6.5, height = 5, dpi = 300)
  p
}

#######################################################################
# boxplot by visit x group, per trait
#######################################################################

plot_trait_boxplot <- function(long_df, trait, label) {
  d <- long_df %>% filter(trait == !!trait) %>% mutate(time = factor(time))
  
  p <- ggplot(d, aes(x = time, y = bmd, fill = Group)) +
    geom_boxplot(outlier.alpha = 0.3, width = 0.6, position = position_dodge(0.7)) +
    scale_fill_manual(values = c(Atenolol = "#D55E00", Placebo = "#0072B2")) +
    labs(title = paste0(label, " by visit and treatment group"),
         x = "Visit month", y = expression(paste("BMD (g/cm"^2, ")"))) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
  
  ggsave(file.path(OUT_DIR, paste0("fig_boxplot_", trait, ".png")),
         plot = p, width = 7, height = 5, dpi = 300)
  p
}

#######################################################################
# MANOVA (basic + full repeated-measures), run per trait
#######################################################################

run_trait_manova <- function(wide_df, trait, times = c(0,6,12,18,24)) {
  cols <- paste0(trait, "_", times)
  missing_cols <- setdiff(cols, names(wide_df))
  if (length(missing_cols) > 0) {
    warning(trait, ": missing time columns ", paste(missing_cols, collapse=", "),
            " -- skipping")
    return(NULL)
  }
  mat <- as.matrix(wide_df[, cols])
  
  cat("\n============================\nMANOVA:", trait, "\n============================\n")
  
  fit <- manova(mat ~ Group * rs1042714, data = wide_df)
  wilks <- summary(fit, test = "Wilks")
  print(wilks)
  
  mod <- lm(mat ~ Group * rs1042714, data = wide_df)
  idata <- data.frame(Time = factor(seq_along(times)))
  rm.manova <- Anova(mod, idata = idata, idesign = ~Time, type = "III")
  full_summary <- summary(rm.manova, multivariate = TRUE)
  print(full_summary)
  
  wilks_stats <- as.data.frame(wilks$stats)
  wilks_stats$term  <- rownames(wilks_stats)
  wilks_stats$trait <- trait
  rownames(wilks_stats) <- NULL
  
  list(basic = fit, full = rm.manova, wilks_table = wilks_stats,
       full_text = capture.output(print(full_summary)))
}

#######################################################################
# SLOPES per subject per trait (mirrors old script's Response/Rebound/Total)
#######################################################################

compute_bmd_slopes <- function(long_df, trait) {
  d <- long_df %>% filter(trait == !!trait, !is.na(time), !is.na(bmd))
  first  <- d %>% filter(time <= 6)
  second <- d %>% filter(time > 6 & time <= 24)
  third  <- d %>% filter(time %in% c(0, 24))
  
  ids <- unique(d$id)
  out <- data.frame(id = ids, Response = NA_real_, Rebound = NA_real_, Total = NA_real_)
  for (i in seq_along(ids)) {
    pid <- ids[i]
    f  <- first[first$id == pid, ]
    s  <- second[second$id == pid, ]
    t3 <- third[third$id == pid, ]
    if (nrow(f)  >= 2) out$Response[i] <- coef(lm(bmd ~ time, data = f))[2]
    if (nrow(s)  >= 2) out$Rebound[i]  <- coef(lm(bmd ~ time, data = s))[2]
    if (nrow(t3) >= 2) out$Total[i]    <- coef(lm(bmd ~ time, data = t3))[2]
  }
  out$trait <- trait
  out
}

#######################################################################
# RUN EVERYTHING, PER TRAIT
#######################################################################

bmd_results <- list()
bmd_slopes_list <- list()

for (trait in unname(BMD_TRAIT_COLS)) {
  label <- BMD_TRAIT_LABELS[[trait]]
  
  # figures
  plot_trait_trajectory(clin_long, trait, label)
  plot_trait_boxplot(clin_long, trait, label)
  
  # MANOVA
  trait_wide <- build_trait_wide(clin_base, trait)
  bmd_results[[trait]] <- run_trait_manova(trait_wide, trait)
  
  # slopes
  bmd_slopes_list[[trait]] <- compute_bmd_slopes(clin_long, trait)
}

bmd_slopes <- bind_rows(bmd_slopes_list) %>%
  left_join(txgroup, by = c("id" = "ID")) %>%
  left_join(geno_wide, by = "id") %>%
  mutate(rs1042714_cat = ifelse(rs1042714 == 0, "GG",
                                ifelse(rs1042714 == 1, "GC", "CC")),
         trait_label = BMD_TRAIT_LABELS[trait])

#######################################################################
# slope boxplots by Group, one figure per slope
# type, faceted across all 5 traits
#######################################################################

plot_slope_boxplot <- function(slopes_df, slope_var, ylab) {
  p <- ggplot(slopes_df, aes(x = Group, y = .data[[slope_var]], fill = Group)) +
    geom_boxplot(outlier.alpha = 0.4, width = 0.6) +
    facet_wrap(~ trait_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Atenolol = "#D55E00", Placebo = "#0072B2")) +
    labs(title = paste0(ylab, " slope by treatment group"), x = NULL,
         y = paste0(ylab, " (g/cm\u00b2 per month)")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey90"),
          axis.text.x = element_text(angle = 0))
  
  ggsave(file.path(OUT_DIR, paste0("fig_slope_", tolower(slope_var), "_by_group.png")),
         plot = p, width = 13, height = 4.5, dpi = 300)
  p
}

plot_slope_boxplot(bmd_slopes, "Response", "Response (0-6mo)")
plot_slope_boxplot(bmd_slopes, "Rebound",  "Rebound (6-24mo)")
plot_slope_boxplot(bmd_slopes, "Total",    "Total (0-24mo)")

#######################################################################
# TABLE: descriptive stats (mean/SD/n) per trait x time x Group
#######################################################################

descriptive_table <- clin_long %>%
  group_by(trait, trait_label, time, Group) %>%
  summarize(n = sum(!is.na(bmd)), mean = mean(bmd, na.rm = TRUE),
            sd = sd(bmd, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean = round(mean, 3), sd = round(sd, 3)) %>%
  arrange(trait, time, Group)

write.csv(descriptive_table, file.path(OUT_DIR, "table1_bmd_descriptive_summary.csv"), row.names = FALSE)

#######################################################################
# TABLE: MANOVA (basic, Wilks) results, all traits,
# rounded + significance stars, Residuals row dropped
#######################################################################

manova_wilks_table <- bind_rows(lapply(bmd_results, function(x) if (!is.null(x)) x$wilks_table)) %>%
  filter(term != "Residuals") %>%
  mutate(trait_label = BMD_TRAIT_LABELS[trait]) %>%
  rename(F_value = `approx F`, p_value = `Pr(>F)`) %>%
  mutate(
    Wilks = round(Wilks, 3), F_value = round(F_value, 2), p_value = round(p_value, 4),
    sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                    p_value < 0.05 ~ "*", TRUE ~ "")
  ) %>%
  select(trait_label, term, Df, Wilks, F_value, `num Df`, `den Df`, p_value, sig)

write.csv(manova_wilks_table, file.path(OUT_DIR, "table2_bmd_manova_wilks.csv"), row.names = FALSE)

full_log <- unlist(lapply(names(bmd_results), function(tr) {
  if (is.null(bmd_results[[tr]])) return(NULL)
  c(paste0("==== ", BMD_TRAIT_LABELS[[tr]], " (", tr, ") ===="),
    bmd_results[[tr]]$full_text, "")
}))
writeLines(full_log, file.path(OUT_DIR, "bmd_manova_full_results.txt"))

#######################################################################
# TABLE: slope Group comparison (t-test), all traits
#######################################################################

slope_types <- c("Response", "Rebound", "Total")
slope_rows <- list()
for (tr in unname(BMD_TRAIT_COLS)) {
  d <- bmd_slopes %>% filter(trait == tr)
  for (sv in slope_types) {
    vals <- d[[sv]]
    grp  <- d$Group
    tt <- tryCatch(t.test(vals ~ grp), error = function(e) NULL)
    slope_rows[[paste(tr, sv)]] <- data.frame(
      trait_label   = BMD_TRAIT_LABELS[[tr]],
      slope_type    = sv,
      n_Atenolol    = sum(grp == "Atenolol" & !is.na(vals)),
      n_Placebo     = sum(grp == "Placebo"  & !is.na(vals)),
      mean_Atenolol = round(mean(vals[grp == "Atenolol"], na.rm = TRUE), 5),
      mean_Placebo  = round(mean(vals[grp == "Placebo"],  na.rm = TRUE), 5),
      p_value       = if (!is.null(tt)) round(tt$p.value, 4) else NA_real_
    )
  }
}
slope_group_table <- bind_rows(slope_rows) %>%
  mutate(sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                         p_value < 0.05 ~ "*", TRUE ~ ""))
write.csv(slope_group_table, file.path(OUT_DIR, "table3_bmd_slope_group_comparison.csv"), row.names = FALSE)

#######################################################################
# TABLE: slope x genotype association (ANOVA), all traits
#######################################################################

geno_rows <- list()
for (tr in unname(BMD_TRAIT_COLS)) {
  d <- bmd_slopes %>% filter(trait == tr)
  for (sv in slope_types) {
    aov_fit <- tryCatch(aov(d[[sv]] ~ d$rs1042714_cat), error = function(e) NULL)
    p <- if (!is.null(aov_fit)) summary(aov_fit)[[1]][["Pr(>F)"]][1] else NA_real_
    geno_rows[[paste(tr, sv)]] <- data.frame(
      trait_label = BMD_TRAIT_LABELS[[tr]], slope_type = sv,
      p_value_rs1042714 = round(p, 4)
    )
  }
}
slope_genotype_table <- bind_rows(geno_rows) %>%
  mutate(sig = case_when(p_value_rs1042714 < 0.001 ~ "***", p_value_rs1042714 < 0.01 ~ "**",
                         p_value_rs1042714 < 0.05 ~ "*", TRUE ~ ""))
write.csv(slope_genotype_table, file.path(OUT_DIR, "table4_bmd_slope_genotype_association.csv"), row.names = FALSE)

