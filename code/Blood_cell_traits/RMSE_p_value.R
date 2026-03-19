library(data.table)

# =========================
# Directories
# =========================
BASE_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"
PHENO_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"
COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates.txt"

traits <- c(
  "WBC_count","RBC_count","Haemoglobin","MCV","RDW",
  "Platelet_count","Plateletcrit","PDW",
  "Lymphocyte_perc","Monocyte_perc",
  "Neutrophill_perc","Eosinophill_perc","Basophill_perc",
  "Reticulocyte_perc","MSCV","HLR_perc"
)

folds <- paste0("fold", 1:5)
pth <- "p10e-8_r0.1"

# =========================
# Load pheno + covar
# =========================
pheno <- fread(PHENO_FN)
pheno[, IID := as.character(IID)]
setkey(pheno, IID)

covar <- fread(COVAR_FN)

# covar 兼容：有 IID 或 FID+IID
if ("IID" %in% names(covar)) {
  covar[, IID := as.character(IID)]
} else if (all(c("FID","IID") %in% names(covar))) {
  covar[, IID := as.character(IID)]
} else {
  stop("COVAR file must contain IID (or FID+IID).")
}
setkey(covar, IID)

# 如需生成 age^2 / age*sex（如果文件里没给）
if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}
if (all(c("age","sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}
if ("sex" %in% names(covar)) covar[, sex := as.factor(sex)]

# =========================
# Choose covariates C
# =========================
# 默认：covar 文件里除 IID/FID 外的所有列都纳入
covar_cols <- setdiff(names(covar), c("FID", "IID"))

# 如果你要严格跟 GWAS 一致，把上面那行换成下面这种（按你真实列名改）
# covar_cols <- c(paste0("pc_genetic", 1:10), "age", "age2", "sex", "age_sex")

# =========================
# Residualize phenotype: y_res = resid(lm(y ~ C))
# =========================
resid_map <- vector("list", length(traits))
names(resid_map) <- traits

for (trait in traits) {
  
  dt_y <- merge(
    pheno[, .(IID, y = get(trait))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  # 去掉缺失（lm 也会处理，但这里先筛掉更明确）
  dt_y <- dt_y[complete.cases(dt_y)]
  if (nrow(dt_y) < 5) {
    warning(sprintf("Too few complete samples for trait %s after merging covariates.", trait))
    next
  }
  
  fml <- as.formula(paste("y ~", paste(covar_cols, collapse = " + ")))
  fit <- lm(fml, data = as.data.frame(dt_y))
  
  resid_map[[trait]] <- data.table(
    IID   = dt_y$IID,
    y_res = as.numeric(residuals(fit))
  )
  setkey(resid_map[[trait]], IID)
}

# =========================
# RMSE function
# =========================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

rmse_list <- list()

for (trait in traits) {
  yres_dt <- resid_map[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) next
  
  for (fold in folds) {
    
    ad_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_beta_clumped.txt")
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn,   col.names = c("IID", "PRS"))
    prs_beta <- fread(beta_fn, col.names = c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    # 用 residual phenotype 合并
    dt_ad <- merge(prs_ad,   yres_dt, by = "IID", all = FALSE)
    dt_be <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
    
    rmse_list[[length(rmse_list) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_ad",
      RMSE  = rmse(dt_ad$y_res, dt_ad$PRS)
    )
    
    rmse_list[[length(rmse_list) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_beta",
      RMSE  = rmse(dt_be$y_res, dt_be$PRS)
    )
  }
}

rmse_dt <- rbindlist(rmse_list)

# =========================
# Mean RMSE (residualized)
# =========================
rmse_mean <- rmse_dt[, .(mean_RMSE = mean(RMSE)), by = .(trait, model)]

rmse_mean_wide <- dcast(rmse_mean, trait ~ model, value.var = "mean_RMSE")
rmse_mean_wide[, better_model := ifelse(PRS_ad < PRS_beta, "PRS_ad", "PRS_beta")]

# =========================
# Paired t-test across folds (residualized RMSE)
# =========================
ttest_res <- rmse_dt[
  , {
    ad <- RMSE[model == "PRS_ad"]
    be <- RMSE[model == "PRS_beta"]
    
    if (length(ad) == length(be) && length(ad) >= 2) {
      tt <- t.test(ad, be, paired = TRUE)
      list(t_stat = unname(tt$statistic), p_value = tt$p.value)
    } else {
      list(t_stat = NA_real_, p_value = NA_real_)
    }
  },
  by = trait
]

final_rmse_summary <- merge(rmse_mean_wide, ttest_res, by = "trait")
final_rmse_summary[, sig := ifelse(p_value < 0.05, "Yes", "No")]
setorder(final_rmse_summary, p_value)

print(final_rmse_summary)

# Significant AND PRS_ad better
sig_ad_better <- final_rmse_summary[p_value < 0.05 & PRS_ad < PRS_beta]
n_sig_ad_better <- nrow(sig_ad_better)

cat("Number of traits where PRS_ad significantly outperforms PRS_beta (p < 0.05) using residualized phenotype:",
    n_sig_ad_better, "\n\n")
print(sig_ad_better)
