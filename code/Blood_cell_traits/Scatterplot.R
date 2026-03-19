library(data.table)
library(ggplot2)

### =========================
### CONFIG
### =========================
BASE_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"
PHENO_FN <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"

traits <- c("Lymphocyte_perc", "RDW", "Platelet_count")
folds <- paste0("fold", 1:5)
pth <- "p10e-8_r0.1"

OUT_DIR <- file.path(BASE_DIR, "Plots", "Scatter_selected_traits_by_model")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

trait_headline <- c(
  "Lymphocyte_perc" = "Lymphocyte_perc: PRS_beta shows stronger predictive performance than PRS_ad",
  "RDW" = "RDW: PRS_ad shows stronger predictive performance than PRS_beta",
  "Platelet_count" = "Platelet_count: No significant difference between PRS_ad and PRS_beta"
)

model_colors <- c("PRS_ad" = "#D73027", "PRS_beta" = "#2C7BB6")

### =========================
### LOAD PHENO
### =========================
pheno <- fread(PHENO_FN)
pheno[, IID := as.character(IID)]

### =========================
### BUILD TRUE VS PREDICTED TABLE
### =========================
scatter_all <- list()

for (trait in traits) {
  for (fold in folds) {
    ad_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_ad.txt")
    beta_fn <- file.path(BASE_DIR, trait, fold, "score", pth, "testing_prs_beta_clumped.txt")
    
    if (!file.exists(ad_fn) || !file.exists(beta_fn)) {
      message("Skipping missing files: ", trait, " / ", fold)
      next
    }
    
    prs_ad <- fread(ad_fn)
    prs_beta <- fread(beta_fn)
    
    setnames(prs_ad, c("IID", "predicted"))
    setnames(prs_beta, c("IID", "predicted"))
    
    prs_ad[, IID := as.character(IID)]
    prs_beta[, IID := as.character(IID)]
    
    true_dt <- pheno[, .(IID, true = get(trait))]
    
    dt_ad <- merge(prs_ad, true_dt, by = "IID")
    dt_beta <- merge(prs_beta, true_dt, by = "IID")
    
    dt_ad[, `:=`(trait = trait, fold = fold, model = "PRS_ad")]
    dt_beta[, `:=`(trait = trait, fold = fold, model = "PRS_beta")]
    
    scatter_all[[length(scatter_all) + 1]] <- dt_ad
    scatter_all[[length(scatter_all) + 1]] <- dt_beta
  }
}

if (length(scatter_all) == 0) {
  stop("No valid trait/fold prediction files were found.")
}

scatter_dt <- rbindlist(scatter_all, use.names = TRUE)
scatter_dt[, trait := factor(trait, levels = traits)]
scatter_dt[, fold := factor(fold, levels = folds)]
scatter_dt[, model := factor(model, levels = c("PRS_ad", "PRS_beta"))]

### =========================
### PLOTS: one figure per trait, faceted by model x fold (2x5 = 10 panels)
### =========================
for (tr in traits) {
  dt_plot <- scatter_dt[trait == tr]
  
  if (nrow(dt_plot) == 0) {
    message("No data for ", tr)
    next
  }
  
  p_trait <- ggplot(dt_plot, aes(x = predicted, y = true, color = model)) +
    geom_point(alpha = 0.28, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.55) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey40"
    ) +
    facet_grid(model ~ fold, scales = "free") +
    scale_color_manual(values = model_colors) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0(trait_headline[[tr]], "\nTrue vs Predicted by PRS model and 5 CV folds"),
      x = "Predicted",
      y = "True",
      color = "PRS type"
    ) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("Scatter_true_vs_pred_", tr, "_model_by_fold_10panels.png")),
    plot = p_trait,
    width = 16,
    height = 7,
    dpi = 300
  )
}

