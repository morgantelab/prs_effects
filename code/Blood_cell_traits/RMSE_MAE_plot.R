library(data.table)
library(ggplot2)

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
# 兼容 covar 文件有 FID/IID 或只有 IID 的情况
if ("IID" %in% names(covar)) {
  covar[, IID := as.character(IID)]
} else if (all(c("FID","IID") %in% names(covar))) {
  covar[, IID := as.character(IID)]
} else {
  stop("COVAR file must contain IID (or FID+IID).")
}
setkey(covar, IID)

# 如果你的 covar 文件里没有这些派生项，就按常见做法生成
# （不会覆盖已存在的列）
if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}
if (all(c("age","sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}

# sex 如果是 0/1/2 这类编码，作为 factor 更符合线性模型里“分类协变量”的常见处理
if ("sex" %in% names(covar)) {
  covar[, sex := as.factor(sex)]
}

### =========================
### METRIC FUNCTIONS
### =========================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2_pred <- function(y, yhat) { cor(y, yhat, use = "complete.obs")^2 }

### =========================
### PRECOMPUTE RESIDUALS: y_res = resid(lm(y ~ C))
### =========================
# 用 covar 文件中的所有协变量列（除 IID/FID）作为 C
covar_cols <- setdiff(names(covar), c("FID", "IID"))

# 如果你想“严格指定”和 GWAS 完全一致的协变量（避免 covar 文件里额外列被纳入），
# 直接在这里手动写死，例如：
# covar_cols <- c(paste0("pc_genetic", 1:10), "age", "age2", "sex", "age_sex")

resid_list <- vector("list", length(traits))
names(resid_list) <- traits

for (trait in traits) {
  dt_y <- merge(
    pheno[, .(IID, y = get(trait))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  # 去掉缺失
  dt_y <- dt_y[complete.cases(dt_y)]
  
  if (nrow(dt_y) == 0) {
    warning(sprintf("No data after merge/complete.cases for trait: %s", trait))
    next
  }
  
  # 拟合 y ~ C
  fml <- as.formula(paste("y ~", paste(covar_cols, collapse = " + ")))
  fit <- lm(fml, data = as.data.frame(dt_y))
  
  resid_list[[trait]] <- data.table(
    IID = dt_y$IID,
    y_res = as.numeric(residuals(fit))
  )
  setkey(resid_list[[trait]], IID)
}

### =========================
### CALC METRICS USING y_res
### =========================
metric_all <- list()

for (trait in traits) {
  yres_dt <- resid_list[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) next
  
  for (fold in folds) {
    
    ad_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_beta_clumped.txt")
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn)
    prs_beta <- fread(beta_fn)
    
    setnames(prs_ad,   c("IID", "PRS"))
    setnames(prs_beta, c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    # merge residual phenotype
    dt_ad <- merge(prs_ad, yres_dt, by = "IID", all = FALSE)
    dt_beta <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
    
    metric_all[[length(metric_all) + 1]] <- data.table(
      trait  = trait,
      fold   = fold,
      model  = "PRS_ad",
      RMSE   = rmse(dt_ad$y_res, dt_ad$PRS),
      MAE    = mae(dt_ad$y_res, dt_ad$PRS),
      R2     = r2_pred(dt_ad$y_res, dt_ad$PRS)
    )
    
    metric_all[[length(metric_all) + 1]] <- data.table(
      trait  = trait,
      fold   = fold,
      model  = "PRS_beta",
      RMSE   = rmse(dt_beta$y_res, dt_beta$PRS),
      MAE    = mae(dt_beta$y_res, dt_beta$PRS),
      R2     = r2_pred(dt_beta$y_res, dt_beta$PRS)
    )
  }
}

metric_dt <- rbindlist(metric_all)

metric_long <- melt(
  metric_dt,
  id.vars = c("trait", "fold", "model"),
  variable.name = "metric",
  value.name = "value"
)

# 固定顺序
metric_long[, trait := factor(trait, levels = traits)]
metric_long[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]
metric_long[, metric := factor(metric, levels = c("RMSE", "MAE"))]

### =========================
### PLOTS (RMSE / MAE / R2)
### =========================
p_rmse <- ggplot(
  metric_long[metric == "RMSE"],
  aes(x = trait, y = value, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  scale_fill_manual(values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV RMSE across 16 traits (residualized phenotype)",
    x = "Trait", y = "RMSE", fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(OUT_DIR, "RMSE_5fold_all_traits_adjusted.png"), p_rmse, width = 10, height = 4, dpi = 300)

p_mae <- ggplot(
  metric_long[metric == "MAE"],
  aes(x = trait, y = value, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  scale_fill_manual(values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV MAE across 16 traits (residualized phenotype)",
    x = "Trait", y = "MAE", fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(OUT_DIR, "MAE_5fold_all_traits_adjusted.png"), p_mae, width = 10, height = 4, dpi = 300)

p_r2 <- ggplot(metric_dt, aes(x = trait, y = R2, fill = model)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  scale_fill_manual(values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV R² across 16 traits (residualized phenotype)",
    x = "Trait", y = "R²", fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(OUT_DIR, "R2_5fold_all_traits_adjusted.png"), p_r2, width = 10, height = 4, dpi = 300)

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
        mean_R2_ad = NA_real_,
        mean_R2_beta = NA_real_,
        diff = NA_real_,
        t_stat = NA_real_,
        p_value = NA_real_
      )
    }
  },
  by = trait
]
r2_ttest[, sig :=
           ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                         ifelse(p_value < 0.05, "*", "ns")))]
r2_ttest[, better :=
           ifelse(mean_R2_ad > mean_R2_beta,
                  "PRS_ad",
                  "PRS_beta")]
print(r2_ttest)

### =========================
### (Optional) Bias plots using residuals: y_res ~ PRS
### =========================
reg_metrics <- list()

for (trait in traits) {
  yres_dt <- resid_list[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) next
  
  for (fold in folds) {
    
    ad_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_beta_clumped.txt")
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn);   setnames(prs_ad,   c("IID", "PRS")); prs_ad[, IID := as.character(IID)]
    prs_beta <- fread(beta_fn); setnames(prs_beta, c("IID", "PRS")); prs_beta[, IID := as.character(IID)]
    
    dt_ad   <- merge(prs_ad,   yres_dt, by = "IID", all = FALSE)
    dt_beta <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
    
    fit_ad   <- lm(y_res ~ PRS, data = dt_ad)
    fit_beta <- lm(y_res ~ PRS, data = dt_beta)
    
    reg_metrics[[length(reg_metrics) + 1]] <- data.table(
      trait = trait, fold = fold, model = "PRS_ad",
      intercept = coef(fit_ad)[1], slope = coef(fit_ad)[2]
    )
    reg_metrics[[length(reg_metrics) + 1]] <- data.table(
      trait = trait, fold = fold, model = "PRS_beta",
      intercept = coef(fit_beta)[1], slope = coef(fit_beta)[2]
    )
  }
}

reg_dt <- rbindlist(reg_metrics)
reg_dt[, trait := factor(trait, levels = traits)]
reg_dt[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]

p_intercept <- ggplot(reg_dt, aes(x = trait, y = intercept, fill = model)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV regression intercept (residual ~ predicted)",
    x = "Trait", y = "Intercept", fill = "PRS type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank())

p_slope <- ggplot(reg_dt, aes(x = trait, y = slope, fill = model)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV regression slope (residual ~ predicted)",
    x = "Trait", y = "Slope", fill = "PRS type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "Intercept_5fold_all_traits_resid.png"), p_intercept, width = 10, height = 4, dpi = 300)
ggsave(file.path(OUT_DIR, "Slope_5fold_all_traits_resid.png"),     p_slope,     width = 10, height = 4, dpi = 300)