#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# 1. Paths
# ============================================================

BLOOD_KVIK_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/kvik"

HEIGHT_BMI_KVIK_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/kvik"

RMSE_FILE <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/poster_18traits_plots/RMSE_paired_ttest_18traits_p10e-8_r0.1.txt"

EFFECT_DIR_BLOOD <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/effects"
EFFECT_DIR_HB    <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/effects"

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/Plots/RMSE_check"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05

# ============================================================
# 2. Trait information
# ============================================================

BLOOD_TRAITS <- c(
  "WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW",
  "Platelet_count", "Plateletcrit", "PDW", "Lymphocyte_perc",
  "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc",
  "Basophill_perc", "Reticulocyte_perc", "MSCV", "HLR_perc"
)

HEIGHT_BMI_TRAITS <- c("HEIGHT", "BMI")

trait_label <- data.table(
  trait = c(
    "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
    "Reticulocyte_perc", "HLR_perc", "Platelet_count", "Plateletcrit",
    "PDW", "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
    "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
    "HEIGHT", "BMI"
  ),
  Trait = c(
    "RBC#", "HGB", "MCV", "RDW", "MSCV",
    "RET%", "HLR%", "PLT#", "PCT",
    "PDW", "WBC#", "LYMPH%", "MONO%",
    "NEUT%", "EO%", "BASO%",
    "HEIGHT", "BMI"
  )
)

trait_order <- trait_label$Trait

# ============================================================
# 3. Shared helper functions
# ============================================================

theme_poster <- function() {
  theme_classic() +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 12, color = "black"),
      legend.text = element_text(size = 11, color = "black"),
      strip.text = element_text(size = 11, face = "bold", color = "black"),
      plot.title = element_blank()
    )
}

make_cor_label <- function(x, y) {
  
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return("Correlation not available")
  }
  
  cor_pearson <- cor.test(x, y, method = "pearson")
  cor_spearman <- cor.test(x, y, method = "spearman", exact = FALSE)
  
  paste0(
    "Pearson r = ", round(unname(cor_pearson$estimate), 3),
    ", P = ", signif(cor_pearson$p.value, 3),
    "\nSpearman rho = ", round(unname(cor_spearman$estimate), 3),
    ", P = ", signif(cor_spearman$p.value, 3)
  )
}

label_position <- function(x, y, x_side = "left", y_side = "top") {
  
  x_rng <- range(x, na.rm = TRUE)
  y_rng <- range(y, na.rm = TRUE)
  
  x_diff <- diff(x_rng)
  y_diff <- diff(y_rng)
  
  if (is.na(x_diff) || x_diff == 0) x_diff <- max(abs(x_rng), na.rm = TRUE)
  if (is.na(y_diff) || y_diff == 0) y_diff <- max(abs(y_rng), na.rm = TRUE)
  if (is.na(x_diff) || x_diff == 0) x_diff <- 1
  if (is.na(y_diff) || y_diff == 0) y_diff <- 1
  
  x_lab <- if (x_side == "left") {
    x_rng[1] + 0.03 * x_diff
  } else {
    x_rng[2] - 0.03 * x_diff
  }
  
  y_lab <- if (y_side == "top") {
    y_rng[2] - 0.04 * y_diff
  } else {
    y_rng[1] + 0.04 * y_diff
  }
  
  list(
    x = x_lab,
    y = y_lab,
    hjust = ifelse(x_side == "left", 0, 1),
    vjust = ifelse(y_side == "top", 1, 0)
  )
}

add_cor_label <- function(label, pos) {
  annotate(
    "label",
    x = pos$x,
    y = pos$y,
    label = label,
    hjust = pos$hjust,
    vjust = pos$vjust,
    size = 4.1,
    label.size = 0.25,
    fill = "white",
    alpha = 0.9
  )
}

detect_col <- function(dt, candidates, label = "column") {
  
  hit <- candidates[candidates %in% names(dt)]
  
  if (length(hit) == 0) {
    stop(
      "Cannot detect ", label, ".\n",
      "Tried: ", paste(candidates, collapse = ", "), "\n",
      "Available columns: ", paste(names(dt), collapse = ", ")
    )
  }
  
  hit[1]
}

save_plot <- function(plot, filename, width = 7.8, height = 5.8) {
  ggsave(
    file.path(OUT_DIR, paste0(filename, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    file.path(OUT_DIR, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height
  )
}

# ============================================================
# 4. Count Bonferroni-significant GWAS-beta / alpha hits
# ============================================================

count_sig_hits <- function(assoc_file, trait, alpha = 0.05) {
  
  if (!file.exists(assoc_file)) {
    warning("Missing assoc file: ", assoc_file)
    return(data.table(
      trait = trait,
      n_gwas_beta_hits = NA_integer_,
      n_tested_snps = NA_integer_,
      bonf_threshold = NA_real_,
      min_p = NA_real_,
      assoc_file = assoc_file
    ))
  }
  
  gwas <- fread(assoc_file)
  
  if (!"Wald_P" %in% names(gwas)) {
    stop(
      "Missing Wald_P column in: ", assoc_file,
      "\nCurrent columns: ", paste(names(gwas), collapse = ", ")
    )
  }
  
  gwas[, Wald_P := as.numeric(Wald_P)]
  gwas <- gwas[!is.na(Wald_P) & Wald_P > 0]
  
  if (nrow(gwas) == 0) {
    warning("No valid Wald_P values in: ", assoc_file)
    return(data.table(
      trait = trait,
      n_gwas_beta_hits = 0L,
      n_tested_snps = 0L,
      bonf_threshold = NA_real_,
      min_p = NA_real_,
      assoc_file = assoc_file
    ))
  }
  
  n_tested <- nrow(gwas)
  bonf_threshold <- alpha / n_tested
  min_p <- min(gwas$Wald_P, na.rm = TRUE)
  n_hits <- sum(gwas$Wald_P < bonf_threshold, na.rm = TRUE)
  
  cat(
    trait,
    "| n_tested =", n_tested,
    "| Bonf =", format(bonf_threshold, scientific = TRUE),
    "| minP =", format(min_p, scientific = TRUE),
    "| hits =", n_hits,
    "\n"
  )
  
  data.table(
    trait = trait,
    n_gwas_beta_hits = n_hits,
    n_tested_snps = n_tested,
    bonf_threshold = bonf_threshold,
    min_p = min_p,
    assoc_file = assoc_file
  )
}

hit_list <- list()

for (i in seq_along(BLOOD_TRAITS)) {
  
  trait <- BLOOD_TRAITS[i]
  
  assoc_file <- file.path(
    BLOOD_KVIK_DIR,
    sprintf("kvik_output.step2.pheno%d.assoc", i)
  )
  
  hit_list[[trait]] <- count_sig_hits(
    assoc_file = assoc_file,
    trait = trait,
    alpha = ALPHA
  )
}

for (trait in HEIGHT_BMI_TRAITS) {
  
  assoc_file <- file.path(
    HEIGHT_BMI_KVIK_DIR,
    trait,
    "kvik_output.step2.assoc"
  )
  
  hit_list[[trait]] <- count_sig_hits(
    assoc_file = assoc_file,
    trait = trait,
    alpha = ALPHA
  )
}

hit_dt <- rbindlist(hit_list, fill = TRUE)

hit_dt <- merge(
  hit_dt,
  trait_label,
  by = "trait",
  all.x = TRUE
)

# ============================================================
# 5. Read RMSE paired t-test result
# ============================================================

rmse_dt <- fread(RMSE_FILE)

required_cols <- c(
  "trait",
  "Trait",
  "mean_RMSE_PGS_ad",
  "mean_RMSE_PGS_beta",
  "mean_diff_ad_minus_beta",
  "significant_bonf"
)

missing_cols <- setdiff(required_cols, names(rmse_dt))

if (length(missing_cols) > 0) {
  stop(
    "RMSE file missing columns: ",
    paste(missing_cols, collapse = ", "),
    "\nCurrent columns: ",
    paste(names(rmse_dt), collapse = ", ")
  )
}

rmse_dt[, mean_RMSE_PGS_ad := as.numeric(mean_RMSE_PGS_ad)]
rmse_dt[, mean_RMSE_PGS_beta := as.numeric(mean_RMSE_PGS_beta)]
rmse_dt[, mean_diff_ad_minus_beta := as.numeric(mean_diff_ad_minus_beta)]

# Positive value means PGS_a,d has lower RMSE
rmse_dt[, RMSE_improvement := mean_RMSE_PGS_beta - mean_RMSE_PGS_ad]

rmse_dt[, check_diff := RMSE_improvement + mean_diff_ad_minus_beta]

if (any(abs(rmse_dt$check_diff) > 1e-8, na.rm = TRUE)) {
  warning("RMSE_improvement does not match -mean_diff_ad_minus_beta. Please check RMSE columns.")
}

# ============================================================
# 6. Merge GWAS-beta hit counts and RMSE improvement
# ============================================================

plot_dt <- merge(
  hit_dt[, .(
    trait,
    Trait_hit = Trait,
    n_gwas_beta_hits,
    n_tested_snps,
    bonf_threshold,
    min_p,
    assoc_file
  )],
  rmse_dt[, .(
    trait,
    Trait,
    mean_RMSE_PGS_ad,
    mean_RMSE_PGS_beta,
    mean_diff_ad_minus_beta,
    RMSE_improvement,
    significant_bonf
  )],
  by = "trait",
  all = FALSE
)

plot_dt <- plot_dt[
  !is.na(n_gwas_beta_hits) &
    !is.na(RMSE_improvement)
]

plot_dt[, Trait := factor(Trait, levels = trait_order)]

cat("\nGWAS-beta hit count + RMSE improvement table:\n")
print(plot_dt[, .(
  trait,
  Trait,
  n_gwas_beta_hits,
  bonf_threshold,
  min_p,
  RMSE_improvement,
  significant_bonf
)])

# ============================================================
# 7. Plot: GWAS-beta hits vs RMSE improvement
# ============================================================

hits_rmse_label <- make_cor_label(
  plot_dt$n_gwas_beta_hits,
  plot_dt$RMSE_improvement
)

hits_rmse_pos <- label_position(
  plot_dt$n_gwas_beta_hits,
  plot_dt$RMSE_improvement,
  x_side = "left",
  y_side = "top"
)

p_hits_rmse <- ggplot(
  plot_dt,
  aes(
    x = n_gwas_beta_hits,
    y = RMSE_improvement
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    aes(shape = significant_bonf),
    size = 3.8,
    alpha = 0.9,
    color = "#4D79A6"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text_repel(
    aes(label = Trait),
    size = 3.6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  add_cor_label(hits_rmse_label, hits_rmse_pos) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("Not significant", "Bonferroni significant")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    x = "Bonferroni-significant GWAS-beta hits",
    y = expression(RMSE[beta] - RMSE[italic(a)*","*italic(d)]),
    shape = "RMSE improvement"
  ) +
  theme_poster()

print(p_hits_rmse)
save_plot(p_hits_rmse, "GWAS_beta_hits_vs_RMSE_improvement", width = 7.8, height = 5.8)

# ============================================================
# 8. Read estimated additive and dominance effects
# ============================================================

read_ad_effects <- function(file_path, trait) {
  
  if (!file.exists(file_path)) {
    warning("Missing a/d effect file: ", file_path)
    return(NULL)
  }
  
  dt <- fread(file_path)
  dt <- as.data.table(dt)
  
  required_base <- c("SNP", "CHR", "POS", "p_d")
  missing_base <- setdiff(required_base, names(dt))
  
  if (length(missing_base) > 0) {
    stop(
      "Missing required columns in ", basename(file_path), ": ",
      paste(missing_base, collapse = ", "),
      "\nCurrent columns: ", paste(names(dt), collapse = ", ")
    )
  }
  
  a_col <- detect_col(
    dt,
    candidates = c(
      "a", "a_hat", "est_a", "estimated_a", "additive", "additive_effect",
      "Additive", "A", "beta_a", "effect_a"
    ),
    label = "estimated additive effect column"
  )
  
  d_col <- detect_col(
    dt,
    candidates = c(
      "d", "d_hat", "est_d", "estimated_d", "dominance", "dominance_effect",
      "Dominance", "D", "beta_d", "effect_d"
    ),
    label = "estimated dominance effect column"
  )
  
  dt[, trait := trait]
  dt[, SNP := as.character(SNP)]
  dt[, CHR := as.integer(as.character(CHR))]
  dt[, POS := as.numeric(as.character(POS))]
  dt[, p_d := as.numeric(p_d)]
  dt[, estimated_a := as.numeric(get(a_col))]
  dt[, estimated_d := as.numeric(get(d_col))]
  
  dt <- dt[
    CHR %in% 1:22 &
      !is.na(SNP) &
      !is.na(CHR) &
      !is.na(POS) &
      !is.na(estimated_a) &
      !is.na(estimated_d)
  ]
  
  if (nrow(dt) == 0) return(NULL)
  
  dt[, .(
    trait,
    SNP,
    CHR,
    POS,
    estimated_a,
    estimated_d,
    p_d
  )]
}

AD_FILE_MAP <- c(
  WBC_count         = file.path(EFFECT_DIR_BLOOD, "WBC_count_ad.txt"),
  RBC_count         = file.path(EFFECT_DIR_BLOOD, "RBC_count_ad.txt"),
  Haemoglobin       = file.path(EFFECT_DIR_BLOOD, "Haemoglobin_ad.txt"),
  MCV               = file.path(EFFECT_DIR_BLOOD, "MCV_ad.txt"),
  RDW               = file.path(EFFECT_DIR_BLOOD, "RDW_ad.txt"),
  Platelet_count    = file.path(EFFECT_DIR_BLOOD, "Platelet_count_ad.txt"),
  Plateletcrit      = file.path(EFFECT_DIR_BLOOD, "Plateletcrit_ad.txt"),
  PDW               = file.path(EFFECT_DIR_BLOOD, "PDW_ad.txt"),
  Lymphocyte_perc   = file.path(EFFECT_DIR_BLOOD, "Lymphocyte_perc_ad.txt"),
  Monocyte_perc     = file.path(EFFECT_DIR_BLOOD, "Monocyte_perc_ad.txt"),
  Neutrophill_perc  = file.path(EFFECT_DIR_BLOOD, "Neutrophill_perc_ad.txt"),
  Eosinophill_perc  = file.path(EFFECT_DIR_BLOOD, "Eosinophill_perc_ad.txt"),
  Basophill_perc    = file.path(EFFECT_DIR_BLOOD, "Basophill_perc_ad.txt"),
  Reticulocyte_perc = file.path(EFFECT_DIR_BLOOD, "Reticulocyte_perc_ad.txt"),
  MSCV              = file.path(EFFECT_DIR_BLOOD, "MSCV_ad.txt"),
  HLR_perc          = file.path(EFFECT_DIR_BLOOD, "HLR_perc_ad.txt"),
  HEIGHT            = file.path(EFFECT_DIR_HB, "HEIGHT_ad.txt"),
  BMI               = file.path(EFFECT_DIR_HB, "BMI_ad.txt")
)

ad_list <- list()

for (trait in names(AD_FILE_MAP)) {
  
  cat("Reading estimated a/d effects:", trait, "\n")
  
  tmp <- read_ad_effects(
    file_path = AD_FILE_MAP[[trait]],
    trait = trait
  )
  
  if (!is.null(tmp) && nrow(tmp) > 0) {
    ad_list[[trait]] <- tmp
  }
}

ad_dt <- rbindlist(ad_list, fill = TRUE)

ad_dt <- merge(
  ad_dt,
  trait_label,
  by = "trait",
  all.x = TRUE
)

ad_dt[, Trait := factor(Trait, levels = trait_order)]

# ============================================================
# 9. Trait-level effect-size summaries
# ============================================================

effect_summary <- ad_dt[
  ,
  .(
    n_effect_snps = .N,
    
    mean_abs_a = mean(abs(estimated_a), na.rm = TRUE),
    median_abs_a = median(abs(estimated_a), na.rm = TRUE),
    mean_a = mean(estimated_a, na.rm = TRUE),
    
    mean_abs_d = mean(abs(estimated_d), na.rm = TRUE),
    median_abs_d = median(abs(estimated_d), na.rm = TRUE),
    mean_d = mean(estimated_d, na.rm = TRUE)
  ),
  by = .(trait, Trait)
]

effect_rmse_dt <- merge(
  effect_summary,
  plot_dt[, .(
    trait,
    RMSE_improvement,
    significant_bonf
  )],
  by = "trait",
  all = FALSE
)

effect_rmse_dt[, Trait := factor(Trait, levels = trait_order)]

cat("\nEffect-size summary merged with RMSE improvement:\n")
print(effect_rmse_dt)

# ============================================================
# 10. Scatterplot distribution of dominance effect size by trait
#     Only plot SNPs with p_d > 0.01
#     x-axis = estimated dominance effect size d
# ============================================================

ad_dt_d_plot <- ad_dt[
  !is.na(p_d) &
    p_d > 0.01 &
    !is.na(estimated_d)
]

p_d_scatter_trait <- ggplot(
  ad_dt_d_plot,
  aes(
    x = estimated_d,
    y = 1
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey50"
  ) +
  geom_jitter(
    width = 0,
    height = 0.18,
    size = 1.8,
    alpha = 0.75,
    color = "#FA7F6F"
  ) +
  facet_wrap(~ Trait, scales = "free_x", ncol = 6) +
  scale_y_continuous(
    breaks = NULL,
    labels = NULL
  ) +
  labs(
    x = "Estimated dominance effect d",
    y = NULL
  ) +
  theme_poster() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 13, color = "black")
  )

print(p_d_scatter_trait)

save_plot(
  p_d_scatter_trait,
  "scatter_distribution_dominance_effect_size_by_trait_pd_gt_0.01",
  width = 12,
  height = 8
)
# ============================================================
# 11. Dominance effect size vs RMSE improvement
# ============================================================

d_rmse_label <- make_cor_label(
  effect_rmse_dt$mean_abs_d,
  effect_rmse_dt$RMSE_improvement
)

d_rmse_pos <- label_position(
  effect_rmse_dt$mean_abs_d,
  effect_rmse_dt$RMSE_improvement,
  x_side = "left",
  y_side = "top"
)

p_d_rmse <- ggplot(
  effect_rmse_dt,
  aes(
    x = mean_abs_d,
    y = RMSE_improvement
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    aes(shape = significant_bonf),
    size = 3.8,
    alpha = 0.9,
    color = "#FA7F6F"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text_repel(
    aes(label = Trait),
    size = 3.6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  add_cor_label(d_rmse_label, d_rmse_pos) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("Not significant", "Bonferroni significant")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    x = "Mean absolute dominance effect size |d|",
    y = expression(RMSE[beta] - RMSE[italic(a)*","*italic(d)]),
    shape = "RMSE improvement"
  ) +
  theme_poster()

print(p_d_rmse)
save_plot(p_d_rmse, "dominance_effect_size_vs_RMSE_improvement", width = 7.8, height = 5.8)

# ============================================================
# 12. Distribution of additive effect size by trait
# ============================================================

p_a_dist_trait <- ggplot(
  ad_dt,
  aes(x = estimated_a)
) +
  geom_histogram(
    bins = 30,
    color = "white",
    fill = "#4D79A6",
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey50"
  ) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 6) +
  labs(
    x = "Estimated additive effect a",
    y = "Number of SNPs"
  ) +
  theme_poster() +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )

print(p_a_dist_trait)
save_plot(p_a_dist_trait, "distribution_additive_effect_size_by_trait", width = 12, height = 8)

# ============================================================
# 13. Additive effect size vs RMSE improvement
# ============================================================

a_rmse_label <- make_cor_label(
  effect_rmse_dt$mean_abs_a,
  effect_rmse_dt$RMSE_improvement
)

a_rmse_pos <- label_position(
  effect_rmse_dt$mean_abs_a,
  effect_rmse_dt$RMSE_improvement,
  x_side = "left",
  y_side = "top"
)

p_a_rmse <- ggplot(
  effect_rmse_dt,
  aes(
    x = mean_abs_a,
    y = RMSE_improvement
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    aes(shape = significant_bonf),
    size = 3.8,
    alpha = 0.9,
    color = "#4D79A6"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text_repel(
    aes(label = Trait),
    size = 3.6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  add_cor_label(a_rmse_label, a_rmse_pos) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("Not significant", "Bonferroni significant")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    x = "Mean absolute additive effect size |a|",
    y = expression(RMSE[beta] - RMSE[italic(a)*","*italic(d)]),
    shape = "RMSE improvement"
  ) +
  theme_poster()

print(p_a_rmse)
save_plot(p_a_rmse, "additive_effect_size_vs_RMSE_improvement", width = 7.8, height = 5.8)

# ============================================================
# 14. SNP-level additive effect size vs dominance effect size
# ============================================================

a_d_label <- make_cor_label(
  ad_dt$estimated_a,
  ad_dt$estimated_d
)

a_d_pos <- label_position(
  ad_dt$estimated_a,
  ad_dt$estimated_d,
  x_side = "left",
  y_side = "top"
)

p_a_d <- ggplot(
  ad_dt,
  aes(
    x = estimated_a,
    y = estimated_d
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey50"
  ) +
  geom_point(
    size = 2.4,
    alpha = 0.8,
    color = "#4D79A6"
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    color = "black"
  ) +
  add_cor_label(a_d_label, a_d_pos) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    x = "Estimated additive effect a",
    y = "Estimated dominance effect d"
  ) +
  theme_poster()

print(p_a_d)
save_plot(p_a_d, "additive_effect_size_vs_dominance_effect_size_SNP_level", width = 7.8, height = 5.8)

cat("\nAll plots saved to:\n")
cat(OUT_DIR, "\n")