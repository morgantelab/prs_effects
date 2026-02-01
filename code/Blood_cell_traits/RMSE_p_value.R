library(data.table)

### =========================
### CONFIG (same as before)
### =========================
BASE_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"
PHENO_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c(
  "WBC_count","RBC_count","Haemoglobin","MCV","RDW",
  "Platelet_count","Plateletcrit","PDW",
  "Lymphocyte_perc","Monocyte_perc",
  "Neutrophill_perc","Eosinophill_perc","Basophill_perc",
  "Reticulocyte_perc","MSCV","HLR_perc"
)

folds <- paste0("fold", 1:5)
pth <- "p10e-8_r0.1"

### =========================
### LOAD PHENO
### =========================
pheno <- fread(PHENO_FN)
pheno[, IID := as.character(IID)]
setkey(pheno, IID)

### =========================
### RMSE FUNCTION
### =========================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

### =========================
### COLLECT RMSE PER FOLD
### =========================
rmse_list <- list()

for (trait in traits) {
  for (fold in folds) {
    
    ad_fn <- file.path(BASE_DIR, trait, fold, "score", pth,
                       "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth,
                         "testing_prs_beta_clumped.txt")
    
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn,   col.names = c("IID", "PRS"))
    prs_beta <- fread(beta_fn, col.names = c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    dt_ad <- merge(prs_ad,   pheno[, .(IID, y = get(trait))], by = "IID")
    dt_be <- merge(prs_beta, pheno[, .(IID, y = get(trait))], by = "IID")
    
    rmse_list[[length(rmse_list) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_ad",
      RMSE  = rmse(dt_ad$y, dt_ad$PRS)
    )
    
    rmse_list[[length(rmse_list) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_beta",
      RMSE  = rmse(dt_be$y, dt_be$PRS)
    )
  }
}

rmse_dt <- rbindlist(rmse_list)

### =========================
### MEAN RMSE (5-fold)
### =========================
rmse_mean <- rmse_dt[
  , .(mean_RMSE = mean(RMSE)), 
  by = .(trait, model)
]

rmse_mean_wide <- dcast(
  rmse_mean,
  trait ~ model,
  value.var = "mean_RMSE"
)

rmse_mean_wide[
  , better_model := ifelse(PRS_ad < PRS_beta, "PRS_ad", "PRS_beta")
]

### =========================
### PAIRED T-TEST (fold-wise)
### =========================
ttest_res <- rmse_dt[
  , {
    ad  <- RMSE[model == "PRS_ad"]
    be  <- RMSE[model == "PRS_beta"]
    
    if (length(ad) == length(be) && length(ad) >= 2) {
      tt <- t.test(ad, be, paired = TRUE)
      list(
        t_stat = unname(tt$statistic),
        p_value = tt$p.value
      )
    } else {
      list(t_stat = NA_real_, p_value = NA_real_)
    }
  },
  by = trait
]

### =========================
### FINAL SUMMARY TABLE
### =========================
final_rmse_summary <- merge(
  rmse_mean_wide,
  ttest_res,
  by = "trait"
)

final_rmse_summary[
  , sig := ifelse(p_value < 0.05, "Yes", "No")
]

setorder(final_rmse_summary, p_value)

print(final_rmse_summary)
