library(data.table)
library(ggplot2)

### =========================
### PATHS
### =========================
PRS_CV_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/prs_cv"

PHENO_FILES <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c("HEIGHT", "BMI")

PLOT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/plots"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

### =========================
### metrics
### =========================
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2_pred <- function(y, yhat) cor(y, yhat, use = "complete.obs")^2

### =========================
### readers (2-column safe)
### =========================
read_prs_2col <- function(fn) {
  x <- fread(fn)
  if (ncol(x) < 2) stop("PRS file has <2 columns: ", fn)
  
  data.table(
    IID = as.character(x[[1]]),
    PRS = as.numeric(x[[2]])
  )
}

list_fold_dirs <- function(trait_dir) {
  if (!dir.exists(trait_dir)) return(character())
  d <- list.dirs(trait_dir, full.names = TRUE, recursive = FALSE)
  d[grepl("^fold\\d+$", basename(d))]
}

list_p_dirs <- function(score_root) {
  if (!dir.exists(score_root)) return(character())
  d <- list.dirs(score_root, full.names = TRUE, recursive = FALSE)
  d[grepl("^p.+_r.+$", basename(d))]
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
  stop("COVAR file must contain IID (or FID+IID): ", COVAR_FN)
}
setkey(covar, IID)

# create age^2 and age*sex if needed
if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}
if (all(c("age", "sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}
if ("sex" %in% names(covar)) {
  covar[, sex := as.factor(sex)]
}

# default: all non-ID columns as covariates
covar_cols <- setdiff(names(covar), c("FID", "IID"))

# If you want to hard-code GWAS covariates, use this instead:
# covar_cols <- c(paste0("pc_genetic", 1:10), "age", "age2", "sex", "age_sex")

### =========================
### PRECOMPUTE RESIDUAL PHENOTYPES
### =========================
resid_map <- vector("list", length(traits))
names(resid_map) <- traits

for (trait in traits) {
  pheno_fn <- PHENO_FILES[[trait]]
  if (!file.exists(pheno_fn)) stop("Missing pheno file: ", pheno_fn)
  
  pheno <- fread(pheno_fn)
  if (!("IID" %in% names(pheno))) stop("Pheno missing IID: ", pheno_fn)
  pheno[, IID := as.character(IID)]
  setkey(pheno, IID)
  
  ph_cols <- setdiff(names(pheno), c("FID", "IID"))
  y_col <- if (trait %in% ph_cols) trait else ph_cols[1]
  if (is.na(y_col) || length(y_col) == 0) {
    stop("No phenotype column in: ", pheno_fn)
  }
  
  dt_y <- merge(
    pheno[, .(IID, y = get(y_col))],
    covar[, c("IID", covar_cols), with = FALSE],
    by = "IID",
    all = FALSE
  )
  
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
  
  cat("[RESID] ", trait, " n=", nrow(resid_map[[trait]]), "\n", sep = "")
}

### =========================
### MAIN: compute per fold, per pdir using residuals
### =========================
metric_all <- list()
reg_all <- list()

for (trait in traits) {
  yres_dt <- resid_map[[trait]]
  if (is.null(yres_dt) || nrow(yres_dt) == 0) {
    cat("[MISS] no residual phenotype for:", trait, "\n")
    next
  }
  
  trait_dir <- file.path(PRS_CV_DIR, trait)
  fold_dirs <- list_fold_dirs(trait_dir)
  if (length(fold_dirs) == 0) {
    cat("[MISS] no fold dirs under:", trait_dir, "\n")
    next
  }
  
  for (fold_dir in fold_dirs) {
    fold <- basename(fold_dir)
    score_root <- file.path(fold_dir, "score")
    p_dirs <- list_p_dirs(score_root)
    
    if (length(p_dirs) == 0) {
      cat("[MISS] no p*_r* under:", score_root, "\n")
      next
    }
    
    for (pdir in p_dirs) {
      pth <- basename(pdir)
      
      ad_fn   <- file.path(pdir, "testing_prs_ad.txt")
      beta_fn <- file.path(pdir, "testing_prs_beta_clumped.txt")
      
      if (!file.exists(ad_fn) || !file.exists(beta_fn)) {
        cat("[MISS] ", trait, " ", fold, " ", pth,
            " ad=", file.exists(ad_fn),
            " beta=", file.exists(beta_fn), "\n", sep = "")
        next
      }
      
      prs_ad   <- read_prs_2col(ad_fn)
      prs_beta <- read_prs_2col(beta_fn)
      
      dt_ad <- merge(prs_ad,   yres_dt, by = "IID", all = FALSE)
      dt_bt <- merge(prs_beta, yres_dt, by = "IID", all = FALSE)
      
      # remove NA just in case
      dt_ad <- dt_ad[!is.na(PRS) & !is.na(y_res)]
      dt_bt <- dt_bt[!is.na(PRS) & !is.na(y_res)]
      
      if (nrow(dt_ad) < 3 || nrow(dt_bt) < 3) {
        cat("[SKIP] too few merged samples: ", trait, " ", fold, " ", pth, "\n", sep = "")
        next
      }
      
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_ad",
        RMSE = rmse(dt_ad$y_res, dt_ad$PRS),
        MAE  = mae(dt_ad$y_res, dt_ad$PRS),
        R2   = r2_pred(dt_ad$y_res, dt_ad$PRS)
      )
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_beta",
        RMSE = rmse(dt_bt$y_res, dt_bt$PRS),
        MAE  = mae(dt_bt$y_res, dt_bt$PRS),
        R2   = r2_pred(dt_bt$y_res, dt_bt$PRS)
      )
      
      fit_ad <- lm(y_res ~ PRS, data = dt_ad)
      fit_bt <- lm(y_res ~ PRS, data = dt_bt)
      
      reg_all[[length(reg_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_ad",
        intercept = unname(coef(fit_ad)[1]),
        slope     = unname(coef(fit_ad)[2])
      )
      reg_all[[length(reg_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_beta",
        intercept = unname(coef(fit_bt)[1]),
        slope     = unname(coef(fit_bt)[2])
      )
      
      cat("[OK] ", trait, " ", fold, " ", pth,
          " n_ad=", nrow(dt_ad), " n_beta=", nrow(dt_bt), "\n", sep = "")
    }
  }
}

metric_dt <- rbindlist(metric_all, fill = TRUE)
reg_dt    <- rbindlist(reg_all, fill = TRUE)

if (nrow(metric_dt) == 0) {
  stop("No metric rows produced. Check [MISS] lines above for missing files/dirs.")
}

metric_long <- melt(
  metric_dt,
  id.vars = c("trait", "fold", "pthresh", "model"),
  variable.name = "metric",
  value.name = "value"
)

metric_long[, trait := factor(trait, levels = traits)]
metric_long[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]
reg_dt[, trait := factor(trait, levels = traits)]

### =========================
### PLOTS
### =========================
plot_box <- function(metric_name, pth, out_png) {
  sub <- metric_long[metric == metric_name & pthresh == pth]
  if (nrow(sub) == 0) return(invisible(NULL))
  
  p <- ggplot(sub, aes(x = trait, y = value, fill = model)) +
    geom_boxplot(
      position = position_dodge(0.8),
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.85,
      linewidth = 0.4
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(metric_name, " (residual phenotype, 5-fold, ", pth, ")"),
      x = "Trait",
      y = metric_name,
      fill = "PRS type"
    ) +
    scale_fill_manual(
      values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(out_png, p, width = 7, height = 4.5, dpi = 300)
  cat("Saved:", out_png, "\n")
}

plot_reg <- function(col, ref_line, pth, out_png) {
  sub <- reg_dt[pthresh == pth]
  if (nrow(sub) == 0) return(invisible(NULL))
  
  p <- ggplot(sub, aes(x = trait, y = get(col), fill = model)) +
    geom_boxplot(
      position = position_dodge(0.8),
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.85,
      linewidth = 0.4
    ) +
    geom_hline(yintercept = ref_line, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(col, " (residual ~ predicted, 5-fold, ", pth, ")"),
      x = "Trait",
      y = col,
      fill = "PRS type"
    ) +
    scale_fill_manual(
      values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(out_png, p, width = 7, height = 4.5, dpi = 300)
  cat("Saved:", out_png, "\n")
}

pths_found <- unique(metric_dt$pthresh)

for (pth in pths_found) {
  plot_box("RMSE", pth, file.path(PLOT_DIR, paste0("RMSE_", pth, "_resid.png")))
  plot_box("MAE",  pth, file.path(PLOT_DIR, paste0("MAE_",  pth, "_resid.png")))
  plot_box("R2",   pth, file.path(PLOT_DIR, paste0("R2_",   pth, "_resid.png")))
  
  plot_reg("intercept", 0, pth, file.path(PLOT_DIR, paste0("Intercept_", pth, "_resid.png")))
  plot_reg("slope",     1, pth, file.path(PLOT_DIR, paste0("Slope_",     pth, "_resid.png")))
}

cat("Done. Saved to: ", PLOT_DIR, "\n", sep = "")

##### T-test
ttest_res <- metric_long[
  ,
  {
    ad  <- value[model == "PRS_ad"]
    be  <- value[model == "PRS_beta"]
    
    if (length(ad) == length(be) && length(ad) >= 2) {
      tt <- t.test(ad, be, paired = TRUE)
      list(
        mean_ad   = mean(ad),
        mean_beta = mean(be),
        diff      = mean(ad - be),
        t_stat    = unname(tt$statistic),
        p_value   = tt$p.value
      )
    } else {
      list(
        mean_ad = NA_real_,
        mean_beta = NA_real_,
        diff = NA_real_,
        t_stat = NA_real_,
        p_value = NA_real_
      )
    }
  },
  by = .(trait, pthresh, metric)
]
ttest_res[, sig :=
            ifelse(p_value < 0.001, "***",
                   ifelse(p_value < 0.01,  "**",
                          ifelse(p_value < 0.05,  "*", "ns")))]
ttest_res[, better :=
            ifelse(metric %in% c("RMSE","MAE"),
                   ifelse(mean_ad < mean_beta, "PRS_ad", "PRS_beta"),
                   ifelse(mean_ad > mean_beta, "PRS_ad", "PRS_beta")
            )]
setorder(ttest_res, metric, p_value)

print(ttest_res)
