library(data.table)
library(ggplot2)

### =========================
### USER CONFIG
### =========================
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"

# IMPORTANT: correlation files are under prs_cv (same as your old script)
BASE_DIR <- file.path(OUT_DIR, "prs_cv")

traits <- c("HEIGHT", "BMI")
folds <- paste0("fold", 1:5)
pthresholds <- c("p10e-8_r0.1", "p10e-5_r0.1")

### =========================
### READ DATA
### =========================
read_one_threshold <- function(pth) {
  res <- list()
  
  for (trait in traits) {
    for (fold in folds) {
      
      dir0 <- file.path(BASE_DIR, trait, fold, "metrics", pth)
      
      beta_fn <- file.path(dir0, "corr_beta.txt")
      ad_fn   <- file.path(dir0, "corr_ad.txt")
      
      if (!file.exists(beta_fn) || !file.exists(ad_fn)) next
      
      beta <- fread(beta_fn, header = FALSE)
      ad   <- fread(ad_fn,   header = FALSE)
      
      res[[length(res) + 1]] <- data.table(
        trait   = trait,
        fold    = fold,
        type    = "PRS_beta",
        corr    = as.numeric(beta[[1]]),
        pthresh = pth
      )
      
      res[[length(res) + 1]] <- data.table(
        trait   = trait,
        fold    = fold,
        type    = "PRS_ad",
        corr    = as.numeric(ad[[1]]),
        pthresh = pth
      )
    }
  }
  
  rbindlist(res, fill = TRUE)
}

dt <- rbindlist(lapply(pthresholds, read_one_threshold), fill = TRUE)

if (nrow(dt) == 0) {
  stop("No correlation files found under: ", BASE_DIR,
       "\nExpected: ", BASE_DIR, "/<trait>/foldX/metrics/<pth>/corr_beta.txt")
}

dt[, trait := factor(trait, levels = traits)]
dt[, type  := factor(type, levels = c("PRS_beta", "PRS_ad"))]

### =========================
### BOXPLOT (2 traits only)
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
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("5-fold CV PRS–phenotype correlation (", pth, ")"),
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
    sd_corr   = sd(corr,   na.rm = TRUE)
  ),
  by = .(trait, type, pthresh)
]

plot_trait_avg_bar <- function(dt_avg, pth) {
  ggplot(
    dt_avg[pthresh == pth],
    aes(x = trait, y = mean_corr, fill = type)
  ) +
    geom_col(position = position_dodge(width = 0.8), width = 0.65, alpha = 0.9) +
    geom_errorbar(
      aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr),
      position = position_dodge(width = 0.8),
      width = 0.25,
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("Average PRS–phenotype correlation (5-fold, ", pth, ")"),
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

# ggsave(file.path(PLOT_DIR, "box_p10e-8_r0.1.png"), p1, width = 7, height = 4.5, dpi = 150)
# ggsave(file.path(PLOT_DIR, "box_p10e-5_r0.1.png"), p2, width = 7, height = 4.5, dpi = 150)
ggsave(file.path(PLOT_DIR, "avgbar_p10e-8_r0.1.png"), p_avg_1, width = 7, height = 4.5, dpi = 150)
ggsave(file.path(PLOT_DIR, "avgbar_p10e-5_r0.1.png"), p_avg_2, width = 7, height = 4.5, dpi = 150)

cat("Done. Saved to: ", PLOT_DIR, "\n", sep = "")
