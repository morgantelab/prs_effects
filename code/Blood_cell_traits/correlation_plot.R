library(data.table)
library(ggplot2)

### =========================
### USER CONFIG
### =========================
BASE_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/prs_cv"

traits <- c(
  "WBC_count",
  "RBC_count",
  "Haemoglobin",
  "MCV",
  "RDW",
  "Platelet_count",
  "Plateletcrit",
  "PDW",
  "Lymphocyte_perc",
  "Monocyte_perc",
  "Neutrophill_perc",
  "Eosinophill_perc",
  "Basophill_perc",
  "Reticulocyte_perc",
  "MSCV",
  "HLR_perc"
)

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
      
      beta <- fread(beta_fn)
      ad   <- fread(ad_fn)
      
      res[[length(res)+1]] <- data.table(
        trait = trait,
        fold  = fold,
        type  = "PRS_beta",
        corr  = beta[[1]],
        pthresh = pth
      )
      
      res[[length(res)+1]] <- data.table(
        trait = trait,
        fold  = fold,
        type  = "PRS_ad",
        corr  = ad[[1]],
        pthresh = pth
      )
    }
  }
  
  rbindlist(res)
}

dt <- rbindlist(lapply(pthresholds, read_one_threshold))

trait_group <- data.table(
  trait = traits,
  panel = rep(c("Traits 1–8", "Traits 9–16"), each = 8)
)

dt <- merge(dt, trait_group, by = "trait", all.x = TRUE)
dt$trait <- factor(dt$trait, levels = traits)
dt$panel <- factor(dt$panel, levels = c("Traits 1–8", "Traits 9–16"))
### =========================
### BOXPLOT FUNCTION
### =========================
plot_trait_box <- function(dt, pth) {
  
  ggplot(
    dt[pthresh == pth],
    aes(x = trait, y = corr, fill = type)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.9),
      width = 0.75,
      linewidth = 0.4,
      outlier.shape = NA,
      alpha = 0.85
    ) +
    scale_fill_manual(
      values = c(
        "PRS_ad"   = "#D73027",
        "PRS_beta" = "#2C7BB6"
      )
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ panel, ncol = 1, scales = "free_x") +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("5-fold CV PRS–phenotype correlation (", pth, ")"),
      x = "Trait",
      y = "Correlation",
      fill = "PRS type"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}



### =========================
### MAKE PLOTS
### =========================
p1 <- plot_trait_box(dt, "p10e-8_r0.1")
p2 <- plot_trait_box(dt, "p10e-5_r0.1")

print(p1)
print(p2)






dt$trait <- factor(dt$trait, levels = traits)

## 5-fold average
dt_avg <- dt[
  , .(
    mean_corr = mean(corr),
    sd_corr   = sd(corr)
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
      width = 0.7,
      alpha = 0.9
    ) +
    geom_errorbar(
      aes(
        ymin = mean_corr - sd_corr,
        ymax = mean_corr + sd_corr
      ),
      position = position_dodge(width = 0.8),
      width = 0.25,
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = c(
        "PRS_ad"   = "#D73027",  # 大红
        "PRS_beta" = "#2C7BB6"   # 蓝
      )
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("Average PRS–phenotype correlation (5-fold CV, ", pth, ")"),
      x = "Trait",
      y = "Mean correlation",
      fill = "PRS type"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p_avg_1 <- plot_trait_avg_bar(dt_avg, "p10e-8_r0.1")
p_avg_2 <- plot_trait_avg_bar(dt_avg, "p10e-5_r0.1")

print(p_avg_1)
print(p_avg_2)


folds <- paste0("fold", 1:5)

pth <- "p10e-8_r0.1"
p_cut <- 0.01

res <- list()

for (trait in traits) {
  for (fold in folds) {
    
    fn <- file.path(
      BASE_DIR,
      trait,
      fold,
      "effects",
      pth,
      "ad_clumped.txt"
    )
    
    if (!file.exists(fn)) {
      res[[length(res) + 1]] <- data.table(
        trait = trait,
        fold  = fold,
        n_dom_snps = NA_integer_
      )
      next
    }
    
    dt <- fread(fn)
    
    ## ====== 如果 p-value 列名不是 p_d，就改这里 ======
    if (!"p_d" %in% names(dt)) {
      stop(paste("Column p_d not found in", fn))
    }
    
    n_dom <- sum(dt$p_d < p_cut, na.rm = TRUE)
    
    res[[length(res) + 1]] <- data.table(
      trait = trait,
      fold  = fold,
      n_dom_snps = n_dom
    )
  }
}

dom_snp_table <- rbindlist(res)
dom_snp_wide <- dcast(
  dom_snp_table,
  trait ~ fold,
  value.var = "n_dom_snps"
)

dom_snp_wide
fold_cols <- paste0("fold", 1:5)

dom_avg <- dom_snp_wide[
  , .(
    mean_dom = rowMeans(.SD, na.rm = TRUE)
  ),
  .SDcols = fold_cols,
  by = trait
]

plot_dt <- merge(
  dt_avg[pthresh == "p10e-8_r0.1"],
  dom_avg,
  by = "trait",
  all.x = TRUE
)
plot_dt$trait <- factor(plot_dt$trait, levels = traits)

corr_max <- max(plot_dt$mean_corr + plot_dt$sd_corr, na.rm = TRUE)
dom_max  <- max(plot_dt$mean_dom, na.rm = TRUE)

scale_factor <- corr_max / dom_max
ggplot(plot_dt, aes(x = trait)) +
  
  ## bars: mean correlation
  geom_col(
    aes(y = mean_corr, fill = type),
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.9
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_corr - sd_corr,
      ymax = mean_corr + sd_corr,
      group = type
    ),
    position = position_dodge(width = 0.8),
    width = 0.25,
    linewidth = 0.4
  ) +
  
  ## line: dominant loci (scaled)
  geom_line(
    aes(y = mean_dom * scale_factor, group = 1),
    color = "black",
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(y = mean_dom * scale_factor),
    color = "black",
    size = 2
  ) +
  
  scale_fill_manual(
    values = c(
      "PRS_ad"   = "#D73027",
      "PRS_beta" = "#2C7BB6"
    )
  ) +
  
  scale_y_continuous(
    name = "Mean PRS–phenotype correlation",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Average number of dominant loci (5-fold)"
    )
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 11) +
  
  labs(
    title = "Average PRS–phenotype correlation with dominant loci usage (p10e-8_r0.1)",
    x = "Trait",
    fill = "PRS type"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
