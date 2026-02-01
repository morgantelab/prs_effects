library(data.table)
library(ggplot2)

## 按你给的顺序固定 trait
traits <- c(
  "WBC_count","RBC_count","Haemoglobin","MCV","RDW",
  "Platelet_count","Plateletcrit","PDW",
  "Lymphocyte_perc","Monocyte_perc",
  "Neutrophill_perc","Eosinophill_perc","Basophill_perc",
  "Reticulocyte_perc","MSCV","HLR_perc"
)

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

  