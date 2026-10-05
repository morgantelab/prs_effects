suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

### =========================
### CONFIG
### =========================

IN_FILE <- paste0(
  "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/",
  "Blood_cell_traits/prs_cv/Plots/",
  "RMSE_R2_5fold_all_traits_adjusted.txt"
)

OUT_DIR <- paste0(
  "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/",
  "Blood_cell_traits/prs_cv/Plots"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


### =========================
### LOAD DATA
### =========================

metric_dt <- fread(IN_FILE)


### =========================
### TRAIT LABEL + ORDER
### =========================

trait_map <- data.table(
  trait = c(
    "RBC_count",
    "Haemoglobin",
    "MCV",
    "RDW",
    "MSCV",
    "Reticulocyte_perc",
    "HLR_perc",
    "Platelet_count",
    "Plateletcrit",
    "PDW",
    "WBC_count",
    "Lymphocyte_perc",
    "Monocyte_perc",
    "Neutrophill_perc",
    "Eosinophill_perc",
    "Basophill_perc"
  ),
  Trait = c(
    "RBC#",
    "HGB",
    "MCV",
    "RDW",
    "MSCV",
    "RET%",
    "HLR%",
    "PLT#",
    "PCT",
    "PDW",
    "WBC#",
    "LYMPH%",
    "MONO%",
    "NEUT%",
    "EO%",
    "BASO%"
  )
)

# Add display labels if not already present
if (!"Trait" %in% names(metric_dt)) {
  metric_dt <- merge(
    metric_dt,
    trait_map,
    by = "trait",
    all.x = TRUE
  )
}

metric_dt[, Trait := factor(
  Trait,
  levels = trait_map$Trait
)]

metric_dt[, model := factor(
  model,
  levels = c("PRS_ad", "PRS_beta")
)]


### =========================
### COLORS + LEGEND LABELS
### =========================

fill_cols <- c(
  "PRS_ad"   = "#FA7F6F",
  "PRS_beta" = "#82B0D2"
)

type_expr <- c(
  expression(plain(PGS)[italic(a) * "," * italic(d)]),
  expression(plain(PGS)[italic(beta)])
)


### =========================
### COMMON THEME
### =========================

my_theme <- theme_classic(base_size = 16) +
  theme(
    text = element_text(
      family = "Arial",
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 14,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 18,
      color = "black",
      margin = margin(r = 8)
    ),
    
    axis.title.x = element_blank(),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.title = element_text(
      size = 16,
      color = "black"
    ),
    
    legend.text = element_text(
      size = 15,
      color = "black"
    ),
    
    legend.key = element_blank(),
    legend.background = element_blank(),
    
    plot.title = element_blank()
  )


### =========================
### R2 PANEL
### =========================

p_r2 <- ggplot(
  metric_dt,
  aes(
    x = Trait,
    y = R2,
    fill = model
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = expression(italic(R)^2)
  ) +
  my_theme +
  theme(
    # Leave blank space on left for manually adding panel label
    plot.margin = margin(
      t = 15,
      r = 10,
      b = 18,
      l = 35
    )
  )


### =========================
### RMSE PANEL
### =========================

p_rmse <- ggplot(
  metric_dt,
  aes(
    x = Trait,
    y = RMSE,
    fill = model
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    linewidth = 0.45,
    outlier.size = 1.3,
    outlier.alpha = 0.6
  ) +
  scale_fill_manual(
    name = "PGS type",
    values = fill_cols,
    breaks = c("PRS_ad", "PRS_beta"),
    labels = type_expr
  ) +
  labs(
    y = "RMSE"
  ) +
  my_theme +
  theme(
    # Leave blank space on left for manually adding panel label
    plot.margin = margin(
      t = 18,
      r = 10,
      b = 15,
      l = 35
    )
  )


### =========================
### COMBINE TWO PANELS
### =========================

combined_plot <- (
  p_r2 /
    p_rmse
) +
  plot_layout(
    ncol = 1,
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )


### =========================
### DISPLAY
### =========================

print(combined_plot)


### =========================
### SAVE
### =========================

ggsave(
  filename = file.path(
    OUT_DIR,
    "R2_RMSE_5fold_combined.png"
  ),
  plot = combined_plot,
  width = 12,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    OUT_DIR,
    "R2_RMSE_5fold_all_traits_adjusted_combined.pdf"
  ),
  plot = combined_plot,
  width = 12,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)
