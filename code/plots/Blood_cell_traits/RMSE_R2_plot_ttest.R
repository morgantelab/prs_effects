suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

### =========================
### CONFIG
### =========================
BASE_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"
PHENO_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"
COVAR    <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates.txt"

traits <- c(
  "WBC_count","RBC_count","Haemoglobin","MCV","RDW",
  "Platelet_count","Plateletcrit","PDW",
  "Lymphocyte_perc","Monocyte_perc",
  "Neutrophill_perc","Eosinophill_perc","Basophill_perc",
  "Reticulocyte_perc","MSCV","HLR_perc"
)

folds <- paste0("fold", 1:5)
pth <- "p10e-8_r0.1"

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv/Plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

### =========================
### LOAD PHENO + COVAR
### =========================
pheno <- fread(PHENO_FN)
pheno[, IID := as.character(IID)]
setkey(pheno, IID)

covar <- fread(COVAR)

if ("IID" %in% names(covar)) {
  covar[, IID := as.character(IID)]
} else if (all(c("FID","IID") %in% names(covar))) {
  covar[, IID := as.character(IID)]
} else {
  stop("COVAR file must contain IID (or FID+IID).")
}
setkey(covar, IID)

if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}
if (all(c("age","sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}
if ("sex" %in% names(covar)) {
  covar[, sex := as.factor(sex)]
}

### =========================
### METRIC FUNCTIONS
### =========================
rmse <- function(y, yhat) {
  sqrt(mean((y - yhat)^2, na.rm = TRUE))
}

r2_pred <- function(y, yhat) {
  cc <- suppressWarnings(cor(y, yhat, use = "complete.obs"))
  if (is.na(cc)) return(NA_real_)
  cc^2
}

### =========================
### PRECOMPUTE RESIDUALS
### =========================
covar_cols <- setdiff(names(covar), c("FID", "IID"))

resid_list <- vector("list", length(traits))
names(resid_list) <- traits

for (trait in traits) {
  dt_y <- merge(
    pheno[, .(IID, y = get(trait))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  dt_y <- dt_y[complete.cases(dt_y)]
  
  if (nrow(dt_y) == 0) {
    warning(sprintf("No data after merge/complete.cases for trait: %s", trait))
    next
  }
  
  fml <- as.formula(paste("y ~", paste(covar_cols, collapse = " + ")))
  fit <- lm(fml, data = as.data.frame(dt_y))
  
  resid_list[[trait]] <- data.table(
    IID   = dt_y$IID,
    y_res = as.numeric(residuals(fit))
  )
  setkey(resid_list[[trait]], IID)
}

### =========================
### CALC METRICS USING RESIDUALS
### =========================
metric_all <- list()

for (trait in traits) {
  yres_dt <- resid_list[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) next
  
  for (fold in folds) {
    ad_fn   <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_beta_clumped.txt")
    
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn)
    prs_beta <- fread(beta_fn)
    
    setnames(prs_ad,   c("IID", "PRS"))
    setnames(prs_beta, c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    prs_ad[, PRS := as.numeric(PRS)]
    prs_beta[, PRS := as.numeric(PRS)]
    
    dt_ad   <- merge(prs_ad,   yres_dt, by = "IID", all = FALSE)
    dt_beta <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
    
    if (nrow(dt_beta) > 1) {
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait  = trait,
        fold   = fold,
        model  = "PRS_beta",
        RMSE   = rmse(dt_beta$y_res, dt_beta$PRS),
        R2     = r2_pred(dt_beta$y_res, dt_beta$PRS)
      )
    }
    
    if (nrow(dt_ad) > 1) {
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait  = trait,
        fold   = fold,
        model  = "PRS_ad",
        RMSE   = rmse(dt_ad$y_res, dt_ad$PRS),
        R2     = r2_pred(dt_ad$y_res, dt_ad$PRS)
      )
    }
  }
}

metric_dt <- rbindlist(metric_all, fill = TRUE)

if (nrow(metric_dt) == 0) {
  stop("No metrics calculated. Check file paths and input files.")
}

fwrite(metric_dt, file.path(OUT_DIR, "RMSE_R2_5fold_all_traits_adjusted.txt"), sep = "\t")

### =========================
### FACTOR ORDER
### =========================
metric_dt[, trait := factor(trait, levels = traits)]
metric_dt[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]

### =========================
### TRAIT LABEL MAP
### =========================
trait_map <- data.table(
  trait_code = c(
    "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
    "Reticulocyte_perc", "HLR_perc", "Platelet_count", "Plateletcrit",
    "PDW", "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
    "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc"
  ),
  Trait = c(
    "RBC#","HGB","MCV","RDW","MSCV","RET%","HLR%",
    "PLT#","PCT","PDW","WBC#","LYMPH%","MONO%",
    "NEUT%","EO%","BASO%"
  )
)

metric_dt <- merge(metric_dt, trait_map, 
                   by.x = "trait", by.y = "trait_code", all.x = TRUE)

metric_dt[, Trait := factor(Trait, levels = trait_map$Trait)]

### =========================
### PLOT STYLE
### =========================
fill_cols <- c(
  "PRS_beta" = "#82B0D2",
  "PRS_ad"   = "#FA7F6F"
)

type_expr <- c(
  expression(plain(PGS)[italic(a) * "," * italic(d)]),
  expression(plain(PGS)[italic(beta)])
)

my_theme <- theme_classic(base_size = 16) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    
    legend.position = "right",
    legend.title = element_text(size = 16, color = "black"),
    legend.text  = element_text(size = 15, color = "black"),
    
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 18,
      color = "black",
      margin = margin(r = 8)
    ),
    axis.title.x = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.key = element_blank(),
    legend.background = element_blank(),
    
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    
    plot.margin = margin(8, 12, 8, 8)
  )
### =========================
### RMSE PLOT: 16 traits in one plot
### =========================
p_rmse <- ggplot(
  metric_dt,
  aes(x = Trait, y = RMSE, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = "RMSE",
    title = NULL
  ) +
  my_theme

p_rmse

ggsave(
  file.path(OUT_DIR, "RMSE_5fold_all_traits_adjusted.png"),
  p_rmse,
  width = 11,
  height = 4.5,
  dpi = 300,
  bg = "white"
)
### =========================
### R2 PLOT: 16 traits in one plot
### =========================
p_r2 <- ggplot(
  metric_dt,
  aes(x = Trait, y = R2, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
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
  my_theme

p_r2
ggsave(
  file.path(OUT_DIR, "R2_5fold_all_traits_adjusted.png"),
  p_r2,
  width = 11,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

### =========================
### PAIRED T-TEST FOR R2
### =========================
r2_ttest <- metric_dt[
  ,
  {
    ad <- R2[model == "PRS_ad"]
    be <- R2[model == "PRS_beta"]
    
    if (length(ad) == length(be) && length(ad) >= 2) {
      tt <- t.test(ad, be, paired = TRUE)
      list(
        mean_R2_ad   = mean(ad),
        mean_R2_beta = mean(be),
        diff         = mean(ad - be),
        t_stat       = unname(tt$statistic),
        p_value      = tt$p.value
      )
    } else {
      list(
        mean_R2_ad   = NA_real_,
        mean_R2_beta = NA_real_,
        diff         = NA_real_,
        t_stat       = NA_real_,
        p_value      = NA_real_
      )
    }
  },
  by = trait
]

r2_ttest[, sig := ifelse(p_value < 0.05/16, "***", "ns")]
r2_ttest[, better := ifelse(mean_R2_ad > mean_R2_beta, "PRS_ad", "PRS_beta")]

print(r2_ttest)
fwrite(r2_ttest, file.path(OUT_DIR, "R2_paired_ttest_5fold_all_traits.txt"), sep = "\t")

### =========================
### PAIRED T-TEST FOR RMSE
### =========================
rmse_ttest <- metric_dt[
  ,
  {
    ad <- RMSE[model == "PRS_ad"]
    be <- RMSE[model == "PRS_beta"]
    
    if (length(ad) == length(be) && length(ad) >= 2) {
      tt <- t.test(ad, be, paired = TRUE)
      list(
        mean_RMSE_ad   = mean(ad),
        mean_RMSE_beta = mean(be),
        diff           = mean(ad - be),
        t_stat         = unname(tt$statistic),
        p_value        = tt$p.value
      )
    } else {
      list(
        mean_RMSE_ad   = NA_real_,
        mean_RMSE_beta = NA_real_,
        diff           = NA_real_,
        t_stat         = NA_real_,
        p_value        = NA_real_
      )
    }
  },
  by = trait
]

rmse_ttest[, sig := ifelse(p_value < 0.05/16, "***", "ns")]
rmse_ttest[, better := ifelse(mean_RMSE_ad < mean_RMSE_beta, "PRS_ad", "PRS_beta")]

print(rmse_ttest)
fwrite(rmse_ttest, file.path(OUT_DIR, "RMSE_paired_ttest_5fold_all_traits.txt"), sep = "\t")
