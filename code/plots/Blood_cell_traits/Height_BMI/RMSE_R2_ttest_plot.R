suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

### =========================
### USER CONFIG
### =========================
OUT_DIR  <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"
BASE_DIR <- file.path(OUT_DIR, "prs_cv")

PHENO_FILES <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c("HEIGHT", "BMI")
folds  <- paste0("fold", 1:5)
pth    <- "p10e-8_r0.1"

PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

### =========================
### HELPER FUNCTIONS
### =========================
calc_rmse <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
}

calc_r2 <- function(y_true, y_pred) {
  cc <- suppressWarnings(cor(y_true, y_pred, use = "complete.obs"))
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
  
  dt <- dt[, .(
    IID = as.character(V1),
    PRS = suppressWarnings(as.numeric(V2))
  )]
  
  dt <- dt[!is.na(IID) & !is.na(PRS)]
  dt
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
  stop("COVAR file must contain IID (or FID+IID).")
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
### LOAD + RESIDUALIZE PHENOTYPES
### =========================
resid_map <- vector("list", length(traits))
names(resid_map) <- traits

for (trait in traits) {
  pheno <- fread(PHENO_FILES[[trait]])
  pheno[, IID := as.character(IID)]
  setkey(pheno, IID)
  
  pheno_cols <- setdiff(names(pheno), c("FID", "IID"))
  if (length(pheno_cols) != 1) {
    stop(sprintf(
      "Phenotype file for %s should contain exactly one phenotype column besides IID/FID. Found: %s",
      trait, paste(pheno_cols, collapse = ", ")
    ))
  }
  
  dt_y <- merge(
    pheno[, .(IID, y = get(pheno_cols))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
  dt_y <- dt_y[complete.cases(dt_y)]
  
  if (nrow(dt_y) < 5) {
    warning(sprintf("Too few complete samples for trait %s", trait))
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

### =========================
### COMPUTE RMSE + R2 BY FOLD
### =========================
metric_all <- list()

for (trait in traits) {
  yres_dt <- resid_map[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) next
  
  for (fold in folds) {
    dir0 <- file.path(BASE_DIR, trait, fold, "score", pth)
    
    beta_fn <- file.path(dir0, "testing_prs_beta_clumped.txt")
    ad_fn   <- file.path(dir0, "testing_prs_ad.txt")
    
    beta <- read_prs(beta_fn)
    ad   <- read_prs(ad_fn)
    
    if (!is.null(beta)) {
      dt_beta <- merge(beta, yres_dt, by = "IID", all = FALSE)
      if (nrow(dt_beta) > 1) {
        metric_all[[length(metric_all) + 1]] <- data.table(
          trait = trait,
          fold  = fold,
          model = "PRS_beta",
          RMSE  = calc_rmse(dt_beta$y_res, dt_beta$PRS),
          R2    = calc_r2(dt_beta$y_res, dt_beta$PRS)
        )
      }
    }
    
    if (!is.null(ad)) {
      dt_ad <- merge(ad, yres_dt, by = "IID", all = FALSE)
      if (nrow(dt_ad) > 1) {
        metric_all[[length(metric_all) + 1]] <- data.table(
          trait = trait,
          fold  = fold,
          model = "PRS_ad",
          RMSE  = calc_rmse(dt_ad$y_res, dt_ad$PRS),
          R2    = calc_r2(dt_ad$y_res, dt_ad$PRS)
        )
      }
    }
  }
}

metric_dt <- rbindlist(metric_all, fill = TRUE)

if (nrow(metric_dt) == 0) {
  stop(
    "No metric results found under: ", BASE_DIR,
    "\nExpected: ", BASE_DIR, "/<trait>/foldX/score/", pth, "/testing_prs_*.txt"
  )
}

fwrite(
  metric_dt,
  file.path(PLOT_DIR, paste0("RMSE_R2_5fold_", pth, "_HEIGHT_BMI.txt")),
  sep = "\t"
)

### =========================
### FACTORS / STYLE
### =========================
trait_label_map <- c(
  HEIGHT = "Height",
  BMI    = "BMI"
)

metric_dt[, trait_label := trait_label_map[trait]]

metric_dt[, trait_label := factor(
  trait_label,
  levels = c("Height", "BMI")
)]

metric_dt[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]

fill_cols <- c(
  "PRS_ad"   = "#FA7F6F",
  "PRS_beta" = "#82B0D2"
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
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black"
    ),
    
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      size = 18,
      color = "black",
      margin = margin(r = 8)
    ),
    
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.key = element_blank(),
    legend.background = element_blank(),
    
    plot.title = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

### =========================
### RMSE PLOT
### =========================
p_rmse <- ggplot(
  metric_dt,
  aes(x = trait_label, y = RMSE, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name   = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = "RMSE"
  ) +
  my_theme

### =========================
### R2 PLOT
### =========================
p_r2 <- ggplot(
  metric_dt,
  aes(x = trait_label, y = R2, fill = model)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name   = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = expression(italic(R)^2)
  ) +
  my_theme

### =========================
### COMBINED R2 + RMSE PLOT
### =========================
p_combined <- (plot_spacer() / p_r2 / plot_spacer() / p_rmse) +
  plot_layout(
    heights = c(0.16, 1, 0.16, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 8)
  )

print(p_combined)

ggsave(
  filename = file.path(PLOT_DIR, paste0("R2_RMSE_combined_boxplot_", pth, "_HEIGHT_BMI.png")),
  plot = p_combined,
  width = 7.5,
  height = 7.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(PLOT_DIR, paste0("R2_RMSE_combined_boxplot_", pth, "_HEIGHT_BMI.pdf")),
  plot = p_combined,
  width = 7.5,
  height = 7.2,
  device = cairo_pdf,
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

r2_ttest[, sig := ifelse(p_value < 0.05/2, "***", "ns")]
r2_ttest[, better := ifelse(mean_R2_ad > mean_R2_beta, "PRS_ad", "PRS_beta")]

print(r2_ttest)

fwrite(
  r2_ttest,
  file.path(PLOT_DIR, paste0("R2_paired_ttest_", pth, "_HEIGHT_BMI.txt")),
  sep = "\t"
)

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

rmse_ttest[, sig := ifelse(p_value < 0.05/2, "***", "ns")]
rmse_ttest[, better := ifelse(mean_RMSE_ad < mean_RMSE_beta, "PRS_ad", "PRS_beta")]

print(rmse_ttest)

fwrite(
  rmse_ttest,
  file.path(PLOT_DIR, paste0("RMSE_paired_ttest_", pth, "_HEIGHT_BMI.txt")),
  sep = "\t"
)

cat("Done. Saved to: ", PLOT_DIR, "\n", sep = "")
