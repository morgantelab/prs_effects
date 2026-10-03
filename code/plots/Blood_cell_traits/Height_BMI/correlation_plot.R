library(data.table)
library(ggplot2)

### =========================
### USER CONFIG
### =========================
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"

BASE_DIR <- file.path(OUT_DIR, "prs_cv")

# separate phenotype files
PHENO_FILES <- list(
  HEIGHT = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  BMI    = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt"
)

# covariate file you provided
COVAR_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c("HEIGHT", "BMI")
folds <- paste0("fold", 1:5)
pthresholds <- c("p10e-8_r0.1", "p10e-5_r0.1")

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

# create derived covariates if needed
if ("age" %in% names(covar) && !"age2" %in% names(covar)) {
  covar[, age2 := as.numeric(age)^2]
}
if (all(c("age", "sex") %in% names(covar)) && !"age_sex" %in% names(covar)) {
  covar[, age_sex := as.numeric(age) * as.numeric(sex)]
}
if ("sex" %in% names(covar)) {
  covar[, sex := as.factor(sex)]
}

# default: use all columns except IID/FID as covariates
covar_cols <- setdiff(names(covar), c("FID", "IID"))

# If you want to hard-code the GWAS covariates instead, use something like:
# covar_cols <- c(paste0("pc_genetic", 1:10), "age", "age2", "sex", "age_sex")

### =========================
### LOAD + RESIDUALIZE PHENOTYPES
### =========================
resid_map <- vector("list", length(traits))
names(resid_map) <- traits

for (trait in traits) {
  pheno <- fread(PHENO_FILES[[trait]])
  pheno[, IID := as.character(IID)]
  setkey(pheno, IID)
  
  # find phenotype column automatically
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
### READ PRS AND RECALCULATE CORR USING y_res
### =========================
read_one_threshold <- function(pth) {
  res <- list()
  
  for (trait in traits) {
    yres_dt <- resid_map[[trait]]
    if (is.null(yres_dt) || nrow(yres_dt) == 0) next
    
    for (fold in folds) {
      dir0 <- file.path(BASE_DIR, trait, fold, "score", pth)
      
      beta_fn <- file.path(dir0, "testing_prs_beta_clumped.txt")
      ad_fn   <- file.path(dir0, "testing_prs_ad.txt")
      
      if (!file.exists(beta_fn) || !file.exists(ad_fn)) next
      
      beta <- fread(beta_fn, header = FALSE, col.names = c("IID", "PRS"))
      ad   <- fread(ad_fn,   header = FALSE, col.names = c("IID", "PRS"))
      
      beta[, IID := as.character(IID)]
      ad[, IID := as.character(IID)]
      
      beta[, PRS := as.numeric(PRS)]
      ad[, PRS := as.numeric(PRS)]
      
      dt_beta <- merge(beta, yres_dt, by = "IID", all = FALSE)
      dt_ad   <- merge(ad,   yres_dt, by = "IID", all = FALSE)
      
      corr_beta <- cor(dt_beta$PRS, dt_beta$y_res, use = "complete.obs")
      corr_ad   <- cor(dt_ad$PRS,   dt_ad$y_res,   use = "complete.obs")
      
      res[[length(res) + 1]] <- data.table(
        trait   = trait,
        fold    = fold,
        type    = "PRS_beta",
        corr    = as.numeric(corr_beta),
        pthresh = pth
      )
      
      res[[length(res) + 1]] <- data.table(
        trait   = trait,
        fold    = fold,
        type    = "PRS_ad",
        corr    = as.numeric(corr_ad),
        pthresh = pth
      )
    }
  }
  
  rbindlist(res, fill = TRUE)
}

dt <- rbindlist(lapply(pthresholds, read_one_threshold), fill = TRUE)

if (nrow(dt) == 0) {
  stop(
    "No PRS files found under: ", BASE_DIR,
    "\nExpected: ", BASE_DIR, "/<trait>/foldX/score/<pth>/testing_prs_*.txt"
  )
}

dt[, trait := factor(trait, levels = traits)]
dt[, type := factor(type, levels = c("PRS_beta", "PRS_ad"))]

### =========================
### BOXPLOT
### =========================
plot_trait_box <- function(dt, pth) {
  ggplot(
    dt[pthresh == pth],
    aes(x = trait, y = corr, fill = type)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.9),
      width = 0.7,
      linewidth = 0.4,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("5-fold CV PRS–residual phenotype correlation (", pth, ")"),
      x = "Trait",
      y = "Correlation",
      fill = "PRS type"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank()
    )
}

p1 <- plot_trait_box(dt, "p10e-8_r0.1")
p2 <- plot_trait_box(dt, "p10e-5_r0.1")

print(p1)
print(p2)

### =========================
### 5-fold average bar
### =========================
dt_avg <- dt[
  , .(
    mean_corr = mean(corr, na.rm = TRUE),
    sd_corr   = sd(corr, na.rm = TRUE)
  ),
  by = .(trait, type, pthresh)
]

plot_trait_avg_bar <- function(dt_avg, pth) {
  ggplot(
    dt_avg[pthresh == pth],
    aes(x = trait, y = mean_corr, fill = type)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.65,
      alpha = 0.9
    ) +
    geom_errorbar(
      aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr),
      position = position_dodge(width = 0.8),
      width = 0.25,
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("Average PRS–residual phenotype correlation (5-fold, ", pth, ")"),
      x = "Trait",
      y = "Mean correlation",
      fill = "PRS type"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank()
    )
}

p_avg_1 <- plot_trait_avg_bar(dt_avg, "p10e-8_r0.1")
p_avg_2 <- plot_trait_avg_bar(dt_avg, "p10e-5_r0.1")

print(p_avg_1)
print(p_avg_2)

### =========================
### SAVE
### =========================
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(PLOT_DIR, "avgbar_p10e-8_r0.1_resid.png"), p_avg_1, width = 7, height = 4.5, dpi = 150)
ggsave(file.path(PLOT_DIR, "avgbar_p10e-5_r0.1_resid.png"), p_avg_2, width = 7, height = 4.5, dpi = 150)

cat("Done. Saved to: ", PLOT_DIR, "\n", sep = "")
