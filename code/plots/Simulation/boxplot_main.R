#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(grid)
})

# ============================================================
# 1. Config
# ============================================================
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal"

KEEP_SCEN <- c("add1000", "mix500_500", "comdom1000")

GRID_ORDER <- c(
  "p10e-3_r0.1", "p10e-3_r0.3", "p10e-3_r0.5",
  "p10e-5_r0.1", "p10e-5_r0.3", "p10e-5_r0.5",
  "p10e-8_r0.1", "p10e-8_r0.3", "p10e-8_r0.5"
)

plot_dir <- file.path(OUT_DIR, "plots_parameter_scan_separate_architectures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Helper functions
# ============================================================
read_pheno <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  dt <- tryCatch(
    fread(path, col.names = c("FID", "IID", "y_true")),
    error = function(e) NULL
  )
  
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  dt[, IID := as.character(IID)]
  dt[, y_true := suppressWarnings(as.numeric(y_true))]
  dt <- dt[!is.na(IID) & !is.na(y_true)]
  
  if (nrow(dt) == 0) return(NULL)
  dt
}

read_prs <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0 || ncol(dt) < 2) return(NULL)
  
  out <- data.table(
    IID    = as.character(dt[[1]]),
    y_pred = suppressWarnings(as.numeric(dt[[2]]))
  )
  
  out <- out[!is.na(IID) & !is.na(y_pred)]
  if (nrow(out) == 0) return(NULL)
  
  out
}

calc_rmse <- function(y_true, y_pred) {
  sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
}

calc_r2 <- function(y_true, y_pred) {
  cc <- suppressWarnings(cor(y_true, y_pred, use = "complete.obs"))
  if (is.na(cc)) return(NA_real_)
  cc^2
}

theme_scan <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = "black"),
      axis.text.x = element_text(
        size = base_size - 2,
        colour = "black",
        angle = 25,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = base_size - 1, colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size - 1, colour = "black"),
      legend.text = element_text(size = base_size - 2, colour = "black"),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5, colour = "black"),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# ============================================================
# 3. Build R2 and RMSE results from raw phenotype + PRS files
# ============================================================
res_list <- list()
idx <- 1L

for (scen in KEEP_SCEN) {
  
  scen_dir <- file.path(OUT_DIR, scen)
  if (!dir.exists(scen_dir)) {
    message("Missing scenario dir: ", scen_dir)
    next
  }
  
  reps <- list.dirs(scen_dir, recursive = FALSE, full.names = FALSE)
  reps <- reps[grepl("^rep", reps)]
  
  if (length(reps) == 0) {
    message("No replicate folders found for: ", scen)
    next
  }
  
  for (rep in reps) {
    
    rep_dir <- file.path(scen_dir, rep)
    
    pheno_file <- file.path(rep_dir, "pheno", "simulated_pheno.txt")
    pheno_dt <- read_pheno(pheno_file)
    
    if (is.null(pheno_dt)) {
      message("Missing or empty phenotype file: ", pheno_file)
      next
    }
    
    for (grid in GRID_ORDER) {
      
      gd <- file.path(rep_dir, "grid", grid)
      if (!dir.exists(gd)) next
      
      beta_file <- file.path(gd, "score", "testing_prs_beta_clumped.txt")
      ad_file   <- file.path(gd, "score", "testing_prs_a_d_method2.txt")
      
      prs_beta <- read_prs(beta_file)
      if (!is.null(prs_beta)) {
        dt_beta <- merge(pheno_dt, prs_beta, by = "IID", all = FALSE)
        
        if (nrow(dt_beta) > 1) {
          res_list[[idx]] <- data.table(
            Scenario = scen,
            Rep      = rep,
            Grid     = grid,
            Type     = "PGS_beta",
            RMSE     = calc_rmse(dt_beta$y_true, dt_beta$y_pred),
            R2       = calc_r2(dt_beta$y_true, dt_beta$y_pred)
          )
          idx <- idx + 1L
        }
      }
      
      prs_ad <- read_prs(ad_file)
      if (!is.null(prs_ad)) {
        dt_ad <- merge(pheno_dt, prs_ad, by = "IID", all = FALSE)
        
        if (nrow(dt_ad) > 1) {
          res_list[[idx]] <- data.table(
            Scenario = scen,
            Rep      = rep,
            Grid     = grid,
            Type     = "PGS_a_d",
            RMSE     = calc_rmse(dt_ad$y_true, dt_ad$y_pred),
            R2       = calc_r2(dt_ad$y_true, dt_ad$y_pred)
          )
          idx <- idx + 1L
        }
      }
    }
  }
}

if (length(res_list) == 0) {
  stop("No metrics were generated. Check file paths and folder structure.")
}

metric_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

# fwrite(
#   metric_dt,
#   file.path(plot_dir, "simulation_R2_RMSE_parameter_scan_source_data.tsv"),
#   sep = "\t"
# )

# ============================================================
# 4. Labels and factor order
# ============================================================
scen_labels <- c(
  add1000    = "Pure additive",
  mix500_500 = "Combined additive and dominance",
  comdom1000 = "Complete dominance"
)

metric_dt <- metric_dt[
  Scenario %in% KEEP_SCEN &
    Grid %in% GRID_ORDER &
    !is.na(RMSE) &
    !is.na(R2)
]

if (nrow(metric_dt) == 0) {
  stop("No data left after filtering.")
}

metric_dt[, Scenario_label := factor(
  Scenario,
  levels = names(scen_labels),
  labels = scen_labels
)]

metric_dt[, Grid := factor(Grid, levels = GRID_ORDER)]
metric_dt[, Type := factor(Type, levels = c("PGS_a_d", "PGS_beta"))]

type_labels <- c(
  PGS_a_d  = expression(plain(PGS)[italic(a) * "," * italic(d)]),
  PGS_beta = expression(plain(PGS)[italic(beta)])
)

fill_cols <- c(
  "PGS_a_d"  = "#FA7F6F",
  "PGS_beta" = "#82B0D2"
)

grid_label_text <- c(
  "p10e-3_r0.1" = "atop(italic(P)*'-value'==10^{-3}, italic(r)^2==0.1)",
  "p10e-3_r0.3" = "atop(italic(P)*'-value'==10^{-3}, italic(r)^2==0.3)",
  "p10e-3_r0.5" = "atop(italic(P)*'-value'==10^{-3}, italic(r)^2==0.5)",
  "p10e-5_r0.1" = "atop(italic(P)*'-value'==10^{-5}, italic(r)^2==0.1)",
  "p10e-5_r0.3" = "atop(italic(P)*'-value'==10^{-5}, italic(r)^2==0.3)",
  "p10e-5_r0.5" = "atop(italic(P)*'-value'==10^{-5}, italic(r)^2==0.5)",
  "p10e-8_r0.1" = "atop(italic(P)*'-value'==10^{-8}, italic(r)^2==0.1)",
  "p10e-8_r0.3" = "atop(italic(P)*'-value'==10^{-8}, italic(r)^2==0.3)",
  "p10e-8_r0.5" = "atop(italic(P)*'-value'==10^{-8}, italic(r)^2==0.5)"
)

# ============================================================
# 5. Plot function: one metric, one architecture, one figure
# ============================================================
make_one_plot <- function(dt, scenario_code, y_col, y_lab, out_name) {
  
  dt_sc <- dt[Scenario == scenario_code]
  if (nrow(dt_sc) == 0) {
    message("No data for scenario: ", scenario_code)
    return(NULL)
  }
  
  p <- ggplot(dt_sc, aes(x = Grid, y = .data[[y_col]], fill = Type)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.62,
      outlier.shape = 16,
      outlier.size = 1.0,
      outlier.alpha = 0.5,
      linewidth = 0.45
    ) +
    labs(
      title = scen_labels[[scenario_code]],
      x = NULL,
      y = y_lab,
      fill = "PGS type"
    ) +
    scale_fill_manual(
      values = fill_cols,
      breaks = c("PGS_a_d", "PGS_beta"),
      labels = type_labels
    ) +
    scale_x_discrete(
      labels = parse(text = grid_label_text[levels(dt_sc$Grid)])
    ) +
    theme_scan(base_size = 14, base_family = "Arial") +
    theme(
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(plot_dir, out_name),
    plot = p,
    width = 11,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  return(p)
}

# ============================================================
# 6. Generate 3 R2 plots + 3 RMSE plots
# ============================================================

# ----- R2 -----
p_r2_add <- make_one_plot(
  metric_dt,
  scenario_code = "add1000",
  y_col = "R2",
  y_lab = expression(italic(R)^2),
  out_name = "simulation_R2_add1000.png"
)

p_r2_mix <- make_one_plot(
  metric_dt,
  scenario_code = "mix500_500",
  y_col = "R2",
  y_lab = expression(italic(R)^2),
  out_name = "simulation_R2_mix500_500.png"
)

p_r2_comdom <- make_one_plot(
  metric_dt,
  scenario_code = "comdom1000",
  y_col = "R2",
  y_lab = expression(italic(R)^2),
  out_name = "simulation_R2_comdom1000.png"
)

# ----- RMSE -----
p_rmse_add <- make_one_plot(
  metric_dt,
  scenario_code = "add1000",
  y_col = "RMSE",
  y_lab = "RMSE",
  out_name = "simulation_RMSE_add1000.png"
)

p_rmse_mix <- make_one_plot(
  metric_dt,
  scenario_code = "mix500_500",
  y_col = "RMSE",
  y_lab = "RMSE",
  out_name = "simulation_RMSE_mix500_500.png"
)

p_rmse_comdom <- make_one_plot(
  metric_dt,
  scenario_code = "comdom1000",
  y_col = "RMSE",
  y_lab = "RMSE",
  out_name = "simulation_RMSE_comdom1000.png"
)

message("Done. Saved 6 figures to: ", plot_dir)
# ============================================================
# 7. Generate one plots
# ============================================================
library(patchwork)

make_one_plot <- function(dt, scenario_code, y_col, y_lab) {
  
  dt_sc <- dt[Scenario == scenario_code]
  
  if (nrow(dt_sc) == 0) {
    message("No data for scenario: ", scenario_code)
    return(NULL)
  }
  
  ggplot(dt_sc, aes(x = Grid, y = .data[[y_col]], fill = Type)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.62,
      outlier.shape = 16,
      outlier.size = 1.0,
      outlier.alpha = 0.5,
      linewidth = 0.45
    ) +
    labs(
      x = NULL,
      y = y_lab,
      fill = "PGS type"
    ) +
    scale_fill_manual(
      values = fill_cols,
      breaks = c("PGS_a_d", "PGS_beta"),
      labels = type_labels
    ) +
    scale_x_discrete(
      labels = parse(text = grid_label_text),
      guide = guide_axis(angle = 45)
    )+
    theme_scan(
      base_size = 14,
      base_family = "Arial"
    ) +
    theme(
      axis.text.x = element_text(
        size = 13,
        angle = 45,
        hjust = 1,
        vjust = 1,
        lineheight = 0.9
      ),
      axis.text.y = element_text(
        size = 13
      ),
      axis.title.y = element_text(
        size = 16
      ),
      legend.title = element_text(
        size = 14
      ),
      legend.text = element_text(
        size = 13
      ),
      plot.tag = element_text(
        face = "bold",
        size = 18
      ),
      plot.tag.position = c(0, 1.03),
      plot.margin = margin(
        t = 25,
        r = 5.5,
        b = 35,
        l = 5.5
      )
    )
}


p_r2_add <- make_one_plot(
  metric_dt, "add1000", "R2", expression(italic(R)^2)
)

p_r2_mix <- make_one_plot(
  metric_dt, "mix500_500", "R2", expression(italic(R)^2)
)

p_r2_comdom <- make_one_plot(
  metric_dt, "comdom1000", "R2", expression(italic(R)^2)
)

panel_tag_theme <- theme(
  plot.tag = element_text(
    face = "bold",
    size = 18,
    hjust = 0,
    vjust = 1
  ),
  
  plot.tag.position = c(0, 1.06),
  
  plot.margin = margin(
    t = 25,
    r = 5.5,
    b = 30,
    l = 5.5
  )
)

p_r2_add <- p_r2_add +
  labs(tag = "A") +
  panel_tag_theme

p_r2_mix <- p_r2_mix +
  labs(tag = "B") +
  panel_tag_theme

p_r2_comdom <- p_r2_comdom +
  labs(tag = "C") +
  panel_tag_theme

p_r2_three <- 
  p_r2_add /
  p_r2_mix /
  p_r2_comdom +
  plot_layout(
    heights = c(
      1,     # panel 1
      1,     # panel 2
      1      # panel 3
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "simulation_R2_three_panel_vertical.png"
  ),
  plot = p_r2_three,
  width = 10,
  height = 15,
  dpi = 300,
  bg = "white"
)

# ----- RMSE -----
p_rmse_add <- make_one_plot(
  metric_dt,
  scenario_code = "add1000",
  y_col = "RMSE",
  y_lab = "RMSE"
)

p_rmse_mix <- make_one_plot(
  metric_dt,
  scenario_code = "mix500_500",
  y_col = "RMSE",
  y_lab = "RMSE"
)

p_rmse_comdom <- make_one_plot(
  metric_dt,
  scenario_code = "comdom1000",
  y_col = "RMSE",
  y_lab = "RMSE"
)

p_rmse_add <- p_rmse_add +
  labs(tag = "A") +
  panel_tag_theme

p_rmse_mix <- p_rmse_mix +
  labs(tag = "B") +
  panel_tag_theme

p_rmse_comdom <- p_rmse_comdom +
  labs(tag = "C") +
  panel_tag_theme

# ----- Combine three RMSE panels vertically -----
p_rmse_three <-
  p_rmse_add /
  p_rmse_mix /
  p_rmse_comdom +
  plot_layout(
    heights = c(
      1,
      1,
      1
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

# ----- Save -----
ggsave(
  filename = file.path(
    plot_dir,
    "simulation_RMSE_three_panel_vertical.png"
  ),
  plot = p_rmse_three,
  width = 10,
  height = 15,
  dpi = 300,
  bg = "white"
)

