suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

### =========================
### CONFIG
### =========================
BASE_DIR_BLOOD <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"
BASE_DIR_HB    <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/prs_cv"

PHENO_BLOOD <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

PHENO_FILES_HB <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates.txt"

folds <- paste0("fold", 1:5)

# only P < 1e-8 and r2 = 0.1
pth <- "p10e-8_r0.1"

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/poster_18traits_plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

### =========================
### TRAITS
### =========================
blood_traits <- c(
  "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
  "Reticulocyte_perc", "HLR_perc",
  "Platelet_count", "Plateletcrit", "PDW",
  "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
  "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc"
)

hb_traits <- c("HEIGHT", "BMI")

all_traits <- c(blood_traits, hb_traits)

trait_map <- data.table(
  trait = c(
    "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
    "Reticulocyte_perc", "HLR_perc",
    "Platelet_count", "Plateletcrit", "PDW",
    "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
    "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
    "HEIGHT", "BMI"
  ),
  Trait = c(
    "RBC#", "HGB", "MCV", "RDW", "MSCV",
    "RET%", "HLR%",
    "PLT#", "PCT", "PDW",
    "WBC#", "LYMPH%", "MONO%",
    "NEUT%", "EO%", "BASO%",
    "HEIGHT", "BMI"
  )
)

### =========================
### FUNCTIONS
### =========================
calc_rmse <- function(y, yhat) {
  sqrt(mean((y - yhat)^2, na.rm = TRUE))
}

calc_r2 <- function(y, yhat) {
  cc <- suppressWarnings(cor(y, yhat, use = "complete.obs"))
  if (is.na(cc)) return(NA_real_)
  cc^2
}

read_prs <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  dt <- tryCatch(
    fread(path, header = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) < 2) return(NULL)
  
  data.table(
    IID = as.character(dt[[1]]),
    PRS = suppressWarnings(as.numeric(dt[[2]]))
  )
}

### =========================
### LOAD COVARIATES
### =========================
covar <- fread(COVAR_FN)

if ("IID" %in% names(covar)) {
  covar[, IID := as.character(IID)]
} else if (all(c("FID", "IID") %in% names(covar))) {
  covar[, IID := as.character(IID)]
} else {
  stop("COVAR file must contain IID or FID+IID.")
}

setkey(covar, IID)

if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}

if (all(c("age", "sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}

if ("sex" %in% names(covar)) {
  covar[, sex := as.factor(sex)]
}

covar_cols <- setdiff(names(covar), c("FID", "IID"))

### =========================
### RESIDUALIZE BLOOD TRAITS
### =========================
resid_map <- list()

pheno_blood <- fread(PHENO_BLOOD)
pheno_blood[, IID := as.character(IID)]
setkey(pheno_blood, IID)

for (trait in blood_traits) {
  if (!trait %in% names(pheno_blood)) {
    warning(sprintf("Trait not found in blood phenotype file: %s", trait))
    next
  }
  
  dt_y <- merge(
    pheno_blood[, .(IID, y = get(trait))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  dt_y <- dt_y[complete.cases(dt_y)]
  
  if (nrow(dt_y) < 5) {
    warning(sprintf("Too few complete samples for trait: %s", trait))
    next
  }
  
  fml <- as.formula(paste("y ~", paste(covar_cols, collapse = " + ")))
  fit <- lm(fml, data = as.data.frame(dt_y))
  
  resid_map[[trait]] <- data.table(
    IID = dt_y$IID,
    y_res = as.numeric(residuals(fit))
  )
  setkey(resid_map[[trait]], IID)
}

### =========================
### RESIDUALIZE HEIGHT + BMI
### =========================
for (trait in hb_traits) {
  pheno <- fread(PHENO_FILES_HB[[trait]])
  pheno[, IID := as.character(IID)]
  setkey(pheno, IID)
  
  pheno_cols <- setdiff(names(pheno), c("FID", "IID"))
  
  if (length(pheno_cols) != 1) {
    stop(sprintf(
      "Phenotype file for %s should contain exactly one phenotype column besides FID/IID. Found: %s",
      trait, paste(pheno_cols, collapse = ", ")
    ))
  }
  
  dt_y <- merge(
    pheno[, .(IID, y = get(pheno_cols[1]))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  dt_y <- dt_y[complete.cases(dt_y)]
  
  if (nrow(dt_y) < 5) {
    warning(sprintf("Too few complete samples for trait: %s", trait))
    next
  }
  
  fml <- as.formula(paste("y ~", paste(covar_cols, collapse = " + ")))
  fit <- lm(fml, data = as.data.frame(dt_y))
  
  resid_map[[trait]] <- data.table(
    IID = dt_y$IID,
    y_res = as.numeric(residuals(fit))
  )
  setkey(resid_map[[trait]], IID)
}

### =========================
### CALCULATE R2 + RMSE
### =========================
metric_all <- list()

for (trait in all_traits) {
  yres_dt <- resid_map[[trait]]
  
  if (is.null(yres_dt) || nrow(yres_dt) == 0) {
    warning(sprintf("No residual phenotype available for trait: %s", trait))
    next
  }
  
  base_dir <- if (trait %in% blood_traits) BASE_DIR_BLOOD else BASE_DIR_HB
  
  for (fold in folds) {
    score_dir <- file.path(base_dir, trait, fold, "score", pth)
    
    ad_fn   <- file.path(score_dir, "testing_prs_ad.txt")
    beta_fn <- file.path(score_dir, "testing_prs_beta_clumped.txt")
    
    prs_ad   <- read_prs(ad_fn)
    prs_beta <- read_prs(beta_fn)
    
    if (!is.null(prs_ad)) {
      dt_ad <- merge(prs_ad, yres_dt, by = "IID", all = FALSE)
      
      if (nrow(dt_ad) > 1) {
        metric_all[[length(metric_all) + 1]] <- data.table(
          trait = trait,
          fold  = fold,
          model = "PRS_ad",
          R2    = calc_r2(dt_ad$y_res, dt_ad$PRS),
          RMSE  = calc_rmse(dt_ad$y_res, dt_ad$PRS)
        )
      }
    }
    
    if (!is.null(prs_beta)) {
      dt_beta <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
      
      if (nrow(dt_beta) > 1) {
        metric_all[[length(metric_all) + 1]] <- data.table(
          trait = trait,
          fold  = fold,
          model = "PRS_beta",
          R2    = calc_r2(dt_beta$y_res, dt_beta$PRS),
          RMSE  = calc_rmse(dt_beta$y_res, dt_beta$PRS)
        )
      }
    }
  }
}

metric_dt <- rbindlist(metric_all, fill = TRUE)

if (nrow(metric_dt) == 0) {
  stop("No metrics calculated. Check PRS file paths and p threshold folder name.")
}

metric_dt <- merge(metric_dt, trait_map, by = "trait", all.x = TRUE)

metric_dt[, Trait := factor(Trait, levels = trait_map$Trait)]
metric_dt[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]

fwrite(
  metric_dt,
  file.path(OUT_DIR, "R2_RMSE_18traits_p10e-8_r0.1.txt"),
  sep = "\t"
)

### =========================
### POSTER STYLE
### =========================
fill_cols <- c(
  "PRS_beta" = "#82B0D2",
  "PRS_ad"   = "#FA7F6F"
)

type_expr <- c(
  expression(plain(PGS)[italic(a) * "," * italic(d)]),
  expression(plain(PGS)[italic(beta)])
)

poster_theme <- theme_classic(base_size = 18) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    
    legend.position = "right",
    legend.title = element_text(size = 18, color = "black"),
    legend.text  = element_text(size = 17, color = "black"),
    
    axis.text.x = element_text(
      size = 17,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 17, color = "black"),
    
    axis.title.y = element_text(size = 21, color = "black"),
    axis.title.x = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5,
      color = "black"
    )
  )

### =========================
### R2 PLOT
### =========================
p_r2 <- ggplot(
  metric_dt,
  aes(x = Trait, y = R2, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.4,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = expression(italic(R)^2),
    title = NULL
  ) +
  poster_theme

print(p_r2)

ggsave(
  filename = file.path(OUT_DIR, "Poster_R2_18traits_p10e-8_r0.1.png"),
  plot = p_r2,
  width = 15,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Poster_R2_18traits_p10e-8_r0.1.pdf"),
  plot = p_r2,
  width = 15,
  height = 6,
  bg = "white"
)

### =========================
### RMSE PLOT
### =========================
### =========================
### RMSE PAIRED T-TEST
### Bonferroni correction: 0.05 / 18
### =========================
alpha_bonf <- 0.05 / 18

rmse_wide <- dcast(
  metric_dt,
  trait + Trait + fold ~ model,
  value.var = "RMSE"
)

rmse_wide <- rmse_wide[
  !is.na(PRS_ad) & !is.na(PRS_beta)
]

rmse_ttest <- rmse_wide[
  ,
  {
    tt <- t.test(PRS_ad, PRS_beta, paired = TRUE)
    
    .(
      n_pairs = .N,
      mean_RMSE_PGS_ad   = mean(PRS_ad, na.rm = TRUE),
      mean_RMSE_PGS_beta = mean(PRS_beta, na.rm = TRUE),
      mean_diff_ad_minus_beta = mean(PRS_ad - PRS_beta, na.rm = TRUE),
      t_statistic = unname(tt$statistic),
      df = unname(tt$parameter),
      p_value = tt$p.value
    )
  },
  by = .(trait, Trait)
]

rmse_ttest[, p_adj_bonf := p.adjust(p_value, method = "bonferroni")]
rmse_ttest[, bonf_threshold := alpha_bonf]
rmse_ttest[, significant_bonf := p_value < alpha_bonf]

setorder(rmse_ttest, Trait)

fwrite(
  rmse_ttest,
  file.path(OUT_DIR, "RMSE_paired_ttest_18traits_p10e-8_r0.1.txt"),
  sep = "\t"
)

print(rmse_ttest)

p_rmse <- ggplot(
  metric_dt,
  aes(x = Trait, y = RMSE, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.4,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = "Prediction Error",
    title = NULL
  ) +
  poster_theme

print(p_rmse)

ggsave(
  filename = file.path(OUT_DIR, "Poster_RMSE_18traits_p10e-8_r0.1.png"),
  plot = p_rmse,
  width = 15,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(OUT_DIR, "Poster_RMSE_18traits_p10e-8_r0.1.pdf"),
  plot = p_rmse,
  width = 15,
  height = 6,
  bg = "white"
)

cat("\nDone.\nSaved files to:\n", OUT_DIR, "\n")