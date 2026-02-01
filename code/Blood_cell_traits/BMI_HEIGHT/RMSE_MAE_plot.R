library(data.table)
library(ggplot2)

### =========================
### PATHS (match your screenshot)
### =========================
PRS_CV_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/prs_cv"

# pheno files are NOT inside PRS_CV_DIR (as you said)
PHENO_FILES <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

traits <- c("HEIGHT", "BMI")

PLOT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/plots"

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
### MAIN: compute per fold, per pdir
### =========================
metric_all <- list()
reg_all <- list()

for (trait in traits) {
  
  pheno_fn <- PHENO_FILES[[trait]]
  if (!file.exists(pheno_fn)) stop("Missing pheno file: ", pheno_fn)
  
  pheno <- fread(pheno_fn)
  if (!("IID" %in% names(pheno))) stop("Pheno missing IID: ", pheno_fn)
  pheno[, IID := as.character(IID)]
  setkey(pheno, IID)
  
  # y column: prefer exact name, else first non-ID column
  ph_cols <- setdiff(names(pheno), c("FID","IID"))
  y_col <- if (trait %in% ph_cols) trait else ph_cols[1]
  if (is.na(y_col) || length(y_col) == 0) stop("No phenotype column in: ", pheno_fn)
  
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
      
      dt_ad <- merge(prs_ad,   pheno[, .(IID, y = get(y_col))], by = "IID")
      dt_bt <- merge(prs_beta, pheno[, .(IID, y = get(y_col))], by = "IID")
      
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_ad",
        RMSE = rmse(dt_ad$y, dt_ad$PRS),
        MAE  = mae(dt_ad$y, dt_ad$PRS),
        R2   = r2_pred(dt_ad$y, dt_ad$PRS)
      )
      metric_all[[length(metric_all) + 1]] <- data.table(
        trait = trait, fold = fold, pthresh = pth, model = "PRS_beta",
        RMSE = rmse(dt_bt$y, dt_bt$PRS),
        MAE  = mae(dt_bt$y, dt_bt$PRS),
        R2   = r2_pred(dt_bt$y, dt_bt$PRS)
      )
      
      fit_ad <- lm(y ~ PRS, data = dt_ad)
      fit_bt <- lm(y ~ PRS, data = dt_bt)
      
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
  id.vars = c("trait","fold","pthresh","model"),
  variable.name = "metric",
  value.name = "value"
)

metric_long[, trait := factor(trait, levels = traits)]
metric_long[, model := factor(model, levels = c("PRS_ad","PRS_beta"))]

plot_box <- function(metric_name, pth, out_png) {
  sub <- metric_long[metric == metric_name & pthresh == pth]
  if (nrow(sub) == 0) return(invisible(NULL))
  
  p <- ggplot(sub, aes(x = trait, y = value, fill = model)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(metric_name, " (5-fold, ", pth, ")"),
      x = "Trait",
      y = metric_name,
      fill = "PRS type"
    ) +
    scale_fill_manual(
      values = c(
        "PRS_ad"   = "#D73027",
        "PRS_beta" = "#2C7BB6"
      )
    ) +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(out_png, p, width = 7, height = 4.5, dpi = 300)
  cat("Saved:", out_png, "\n")
}
reg_dt[, trait := factor(trait, levels = c("HEIGHT", "BMI"))]


plot_reg <- function(col, ref_line, pth, out_png) {
  sub <- reg_dt[pthresh == pth]
  if (nrow(sub) == 0) return(invisible(NULL))
  
  p <- ggplot(sub, aes(x = trait, y = get(col), fill = model)) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.85, linewidth = 0.4) +
    geom_hline(yintercept = ref_line, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(col, " (true ~ predicted, 5-fold, ", pth, ")"),
      x = "Trait",
      y = col,
      fill = "PRS type"
    ) +
    scale_fill_manual(
      values = c(
        "PRS_ad"   = "#D73027",
        "PRS_beta" = "#2C7BB6"
      )
    ) +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(out_png, p, width = 7, height = 4.5, dpi = 300)
  cat("Saved:", out_png, "\n")
}

pths_found <- unique(metric_dt$pthresh)

for (pth in pths_found) {
  plot_box("RMSE", pth, file.path(PLOT_DIR, paste0("RMSE_", pth, ".png")))
  plot_box("MAE",  pth, file.path(PLOT_DIR, paste0("MAE_",  pth, ".png")))
  plot_box("R2",   pth, file.path(PLOT_DIR, paste0("R2_",   pth, ".png")))
  
  plot_reg("intercept", 0, pth, file.path(PLOT_DIR, paste0("Intercept_", pth, ".png")))
  plot_reg("slope",     1, pth, file.path(PLOT_DIR, paste0("Slope_",     pth, ".png")))
}

cat("Done. Saved to: ", PLOT_DIR, "\n", sep = "")
