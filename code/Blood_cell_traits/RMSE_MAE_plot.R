library(data.table)
library(ggplot2)

### =========================
### CONFIG
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

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv/Plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

### =========================
### LOAD PHENO
### =========================
pheno <- fread(PHENO_FN)
pheno$IID <- as.character(pheno$IID)
setkey(pheno, IID)

### =========================
### METRIC FUNCTIONS
### =========================

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2_pred <- function(y, yhat) {cor(y, yhat, use = "complete.obs")^2}

metric_all <- list()

for (trait in traits) {
  for (fold in folds) {
    
    ad_fn <- file.path(
      BASE_DIR, trait, fold, "score", pth,
      "testing_prs_ad.txt"
    )
    
    beta_fn <- file.path(
      BASE_DIR, trait, fold, "score", pth,
      "testing_prs_beta_clumped.txt"
    )
    
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn)
    prs_beta <- fread(beta_fn)
    
    setnames(prs_ad,   c("IID", "PRS"))
    setnames(prs_beta, c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    dt_ad <- merge(
      prs_ad,
      pheno[, .(IID, y = get(trait))],
      by = "IID"
    )
    
    dt_beta <- merge(
      prs_beta,
      pheno[, .(IID, y = get(trait))],
      by = "IID"
    )
    
    metric_all[[length(metric_all) + 1]] <- data.table(
      trait  = trait,
      fold   = fold,
      model  = "PRS_ad",
      RMSE   = rmse(dt_ad$y, dt_ad$PRS),
      MAE    = mae(dt_ad$y, dt_ad$PRS),
      R2     = r2_pred(dt_ad$y, dt_ad$PRS)
    )
    
    metric_all[[length(metric_all) + 1]] <- data.table(
      trait  = trait,
      fold   = fold,
      model  = "PRS_beta",
      RMSE   = rmse(dt_beta$y, dt_beta$PRS),
      MAE    = mae(dt_beta$y, dt_beta$PRS),
      R2     = r2_pred(dt_beta$y, dt_beta$PRS)
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

## 固定顺序
metric_long$trait  <- factor(metric_long$trait, levels = traits)
metric_long$model  <- factor(metric_long$model, levels = c("PRS_ad", "PRS_beta"))
metric_long$metric <- factor(metric_long$metric, levels = c("RMSE", "MAE"))

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
  scale_fill_manual(
    values = c(
      "PRS_ad"   = "#D73027",
      "PRS_beta" = "#2C7BB6"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV RMSE across 16 traits",
    x = "Trait",
    y = "RMSE",
    fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "RMSE_5fold_all_traits.png"),
  plot = p_rmse,
  width = 10,
  height = 4,
  dpi = 300
)


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
  scale_fill_manual(
    values = c(
      "PRS_ad"   = "#D73027",
      "PRS_beta" = "#2C7BB6"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV MAE across 16 traits",
    x = "Trait",
    y = "MAE",
    fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "MAE_5fold_all_traits.png"),
  plot = p_mae,
  width = 10,
  height = 4,
  dpi = 300
)

### Bias plot
reg_metrics <- list()

for (trait in traits) {
  for (fold in folds) {
    
    ad_fn <- file.path(
      BASE_DIR, trait, fold, "score", pth,
      "testing_prs_ad.txt"
    )
    
    beta_fn <- file.path(
      BASE_DIR, trait, fold, "score", pth,
      "testing_prs_beta_clumped.txt"
    )
    
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) next
    
    prs_ad   <- fread(ad_fn)
    prs_beta <- fread(beta_fn)
    
    setnames(prs_ad,   c("IID", "PRS"))
    setnames(prs_beta, c("IID", "PRS"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    dt_ad <- merge(
      prs_ad,
      pheno[, .(IID, y = get(trait))],
      by = "IID"
    )
    
    dt_beta <- merge(
      prs_beta,
      pheno[, .(IID, y = get(trait))],
      by = "IID"
    )
    
    ## linear regression: true ~ predicted
    fit_ad   <- lm(y ~ PRS, data = dt_ad)
    fit_beta <- lm(y ~ PRS, data = dt_beta)
    
    reg_metrics[[length(reg_metrics) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_ad",
      intercept = coef(fit_ad)[1],
      slope     = coef(fit_ad)[2]
    )
    
    reg_metrics[[length(reg_metrics) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      model = "PRS_beta",
      intercept = coef(fit_beta)[1],
      slope     = coef(fit_beta)[2]
    )
  }
}

reg_dt<- rbindlist(reg_metrics)

reg_dt$trait <- factor(reg_dt$trait, levels = traits)
reg_dt$model <- factor(reg_dt$model, levels = c("PRS_ad", "PRS_beta"))

p_intercept <- ggplot(
  reg_dt,
  aes(x = trait, y = intercept, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "PRS_ad"   = "#D73027",
      "PRS_beta" = "#2C7BB6"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV regression intercept (true ~ predicted)",
    x = "Trait",
    y = "Intercept",
    fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
  )

p_slope <- ggplot(
  reg_dt,
  aes(x = trait, y = slope, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "PRS_ad"   = "#D73027",
      "PRS_beta" = "#2C7BB6"
    )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV regression slope (true ~ predicted)",
    x = "Trait",
    y = "Slope",
    fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(OUT_DIR, "Intercept_5fold_all_traits.png"),
  plot = p_intercept,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(OUT_DIR, "Slope_5fold_all_traits.png"),
  plot = p_slope,
  width = 10,
  height = 4,
  dpi = 300
)


p_r2 <- ggplot(metric_dt, aes(x = trait, y = R2, fill = model)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "5-fold CV R² across 16 traits",
    x = "Trait",
    y = "R²",
    fill = "PRS type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "R2_5fold_all_traits.png"),
  plot = p_r2,
  width = 10,
  height = 4,
  dpi = 300
)
