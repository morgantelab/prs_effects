suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

### =========================
### USER CONFIG
### =========================
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"
BASE_DIR <- file.path(OUT_DIR, "prs_cv")

PHENO_FILES <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c("HEIGHT", "BMI")
folds <- paste0("fold", 1:5)
pthresholds <- c("p10e-8_r0.1", "p10e-5_r0.1")

plot_dir <- file.path(OUT_DIR, "accuracy_plots_residualized")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

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
read_one_threshold_metrics <- function(pth) {
  res <- list()
  
  for (trait in traits) {
    yres_dt <- resid_map[[trait]]
    if (is.null(yres_dt) || nrow(yres_dt) == 0) next
    
    for (fold in folds) {
      dir0 <- file.path(BASE_DIR, trait, fold, "score", pth)
      
      beta_fn <- file.path(dir0, "testing_prs_beta_clumped.txt")
      ad_fn   <- file.path(dir0, "testing_prs_ad.txt")
      
      # beta
      if (file.exists(beta_fn)) {
        beta <- read_prs(beta_fn)
        if (!is.null(beta)) {
          dt_beta <- merge(beta, yres_dt, by = "IID", all = FALSE)
          if (nrow(dt_beta) > 1) {
            rmse_beta <- calc_rmse(dt_beta$y_res, dt_beta$PRS)
            r2_beta   <- calc_r2(dt_beta$y_res, dt_beta$PRS)
            
            res[[length(res) + 1]] <- data.table(
              trait   = trait,
              fold    = fold,
              type    = "PRS_beta",
              rmse    = as.numeric(rmse_beta),
              r2      = as.numeric(r2_beta),
              pthresh = pth
            )
          }
        }
      }
      
      # ad
      if (file.exists(ad_fn)) {
        ad <- read_prs(ad_fn)
        if (!is.null(ad)) {
          dt_ad <- merge(ad, yres_dt, by = "IID", all = FALSE)
          if (nrow(dt_ad) > 1) {
            rmse_ad <- calc_rmse(dt_ad$y_res, dt_ad$PRS)
            r2_ad   <- calc_r2(dt_ad$y_res, dt_ad$PRS)
            
            res[[length(res) + 1]] <- data.table(
              trait   = trait,
              fold    = fold,
              type    = "PRS_ad",
              rmse    = as.numeric(rmse_ad),
              r2      = as.numeric(r2_ad),
              pthresh = pth
            )
          }
        }
      }
    }
  }
  
  rbindlist(res, fill = TRUE)
}

metric_dt <- rbindlist(lapply(pthresholds, read_one_threshold_metrics), fill = TRUE)

if (nrow(metric_dt) == 0) {
  stop(
    "No metric results found under: ", BASE_DIR,
    "\nExpected: ", BASE_DIR, "/<trait>/foldX/score/<pth>/testing_prs_*.txt"
  )
}

fwrite(metric_dt, file.path(plot_dir, "RMSE_R2_by_fold_residualized.txt"), sep = "\t")

### =========================
### FACTORS / LABELS
### =========================
metric_dt[, trait := factor(trait, levels = traits)]
metric_dt[, type := factor(type, levels = c("PRS_ad", "PRS_beta"))]

fill_cols <- c(
  "PRS_ad"   = "#FA7F6F",
  "PRS_beta" = "#82B0D2"
)

type_expr <- c(
  expression(PGS[a * "," * d]),
  expression(PGS[beta])
)

### =========================
### COMMON THEME
### =========================
my_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11),
    
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    
    strip.text = element_text(face = "bold", size = 12),
    
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
  )

### =========================
### PLOT FUNCTIONS
### =========================
plot_trait_rmse <- function(dt, pth) {
  ggplot(
    dt[pthresh == pth],
    aes(x = type, y = rmse, fill = type)
  ) +
    geom_boxplot(
      width = 0.65,
      linewidth = 0.4,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    facet_wrap(~ trait, nrow = 1, scales = "free_y") +
    scale_fill_manual(
      values = fill_cols,
      breaks = c("PRS_ad", "PRS_beta"),
      labels = type_expr
    ) +
    scale_x_discrete(
      breaks = c("PRS_ad", "PRS_beta"),
      labels = type_expr
    ) +
    labs(
      x = NULL,
      y = "RMSE",
      fill = NULL
    ) +
    my_theme
}

plot_trait_r2 <- function(dt, pth) {
  ggplot(
    dt[pthresh == pth],
    aes(x = type, y = r2, fill = type)
  ) +
    geom_boxplot(
      width = 0.65,
      linewidth = 0.4,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    facet_wrap(~ trait, nrow = 1, scales = "free_y") +
    scale_fill_manual(
      values = fill_cols,
      breaks = c("PRS_ad", "PRS_beta"),
      labels = type_expr
    ) +
    scale_x_discrete(
      breaks = c("PRS_ad", "PRS_beta"),
      labels = type_expr
    ) +
    labs(
      x = NULL,
      y = expression(R^2),
      fill = NULL
    ) +
    my_theme
}

### =========================
### MAKE PLOTS
### =========================
for (pth in pthresholds) {
  p_rmse <- plot_trait_rmse(metric_dt, pth)
  p_r2   <- plot_trait_r2(metric_dt, pth)
  
  print(p_rmse)
  print(p_r2)
  
  ggsave(
    filename = file.path(plot_dir, paste0("RMSE_boxplot_", pth, ".png")),
    plot = p_rmse,
    width = 8.5,
    height = 4.8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0("R2_boxplot_", pth, ".png")),
    plot = p_r2,
    width = 8.5,
    height = 4.8,
    dpi = 300
  )
}
