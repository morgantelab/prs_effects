suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(grid)
})

# ============================================================
# 1. User config
# ============================================================
base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal"

KEEP_SCEN <- c("add1000", "mix500_500", "comdom1000")

plot_dir <- file.path(base_dir, "plots_supplementary_R2_RMSE")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Helpers
# ============================================================
read_one_value <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NA_real_)
  suppressWarnings(as.numeric(dt[[1]][1]))
}

theme_supp <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x   = element_text(
        size = base_size - 2,
        colour = "black",
        angle = 20,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y   = element_text(size = base_size - 1, colour = "black"),
      axis.title.x  = element_text(size = base_size, colour = "black"),
      axis.title.y  = element_text(size = base_size, colour = "black"),
      legend.title  = element_text(size = base_size - 1, colour = "black"),
      legend.text   = element_text(size = base_size - 2, colour = "black"),
      axis.line     = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks    = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.key    = element_blank(),
      legend.background = element_blank(),
      panel.border  = element_blank(),
      plot.margin   = margin(8, 10, 8, 8)
    )
}

# ============================================================
# 3. Labels and factor order
# ============================================================
scen_labels <- c(
  add1000    = "Pure additive",
  mix500_500 = "Combined additive and dominance",
  comdom1000 = "Complete dominance"
)

grid_order <- c(
  "p10e-3_r0.1", "p10e-3_r0.3", "p10e-3_r0.5",
  "p10e-5_r0.1", "p10e-5_r0.3", "p10e-5_r0.5",
  "p10e-8_r0.1", "p10e-8_r0.3", "p10e-8_r0.5"
)

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
# 4. Function to read metric files
# ============================================================
read_metric_files <- function(base_dir, metric = c("R2", "RMSE")) {
  
  metric <- match.arg(metric)
  
  if (metric == "R2") {
    file_pattern <- "^corr_(beta|ad)\\.txt$"
  } else {
    file_pattern <- "^rmse_(beta|ad)\\.txt$"
  }
  
  files <- list.files(
    base_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop("No ", metric, " files found. Check base_dir, folder structure, or file names.")
  }
  
  res <- rbindlist(lapply(files, function(f) {
    
    parts <- str_split(f, "/", simplify = TRUE)
    
    # Folder structure assumed from your original code:
    # ... / Scenario / Rep / ... / Grid / ... / metric_file
    scen <- parts[ncol(parts) - 5]
    rep  <- parts[ncol(parts) - 4]
    grid <- parts[ncol(parts) - 2]
    
    pval <- str_match(grid, "p([^_]+)")[, 2]
    ldv  <- str_match(grid, "r(.+)")[, 2]
    
    if (metric == "R2") {
      type <- ifelse(str_detect(basename(f), "corr_beta"), "PGS_beta", "PGS_a_d")
      val_raw <- read_one_value(f)
      val <- val_raw^2
    } else {
      type <- ifelse(str_detect(basename(f), "rmse_beta"), "PGS_beta", "PGS_a_d")
      val <- read_one_value(f)
    }
    
    data.table(
      Scenario = scen,
      Rep      = rep,
      P        = paste0("10E-", sub(".*e-", "", pval)),
      LD       = ldv,
      Grid     = grid,
      Type     = type,
      Metric   = metric,
      Value    = val
    )
  }), fill = TRUE)
  
  res <- res[
    Scenario %in% KEEP_SCEN & !is.na(Value)
  ]
  
  res <- res[Grid %in% grid_order]
  
  if (nrow(res) == 0) {
    stop("No data remain for ", metric, " after filtering. Check Scenario names and Grid names.")
  }
  
  res[, Scenario_label := factor(
    Scenario,
    levels = names(scen_labels),
    labels = scen_labels
  )]
  
  res[, Type := factor(Type, levels = c("PGS_a_d", "PGS_beta"))]
  res[, Grid := factor(Grid, levels = grid_order)]
  
  return(res)
}

# ============================================================
# 5. Function to plot one metric
# ============================================================
plot_metric_by_architecture <- function(dt, metric = c("R2", "RMSE")) {
  
  metric <- match.arg(metric)
  
  if (metric == "R2") {
    y_lab <- expression(italic(R)^2)
    file_prefix <- "PGS_R2_supplementary"
  } else {
    y_lab <- "RMSE"
    file_prefix <- "PGS_RMSE_supplementary"
  }
  
  for (sc in KEEP_SCEN) {
    
    dt_sc <- dt[Scenario == sc]
    if (nrow(dt_sc) == 0) next
    
    p_sc <- ggplot(dt_sc, aes(x = Grid, y = Value, fill = Type)) +
      geom_boxplot(
        position = position_dodge(width = 0.75),
        width = 0.62,
        outlier.shape = 16,
        outlier.size = 1.0,
        outlier.alpha = 0.5,
        linewidth = 0.45
      ) +
      labs(
        title = scen_labels[sc],
        x = NULL,
        y = y_lab,
        fill = "PGS type"
      ) +
      scale_fill_manual(
        values = fill_cols,
        labels = type_labels
      ) +
      scale_x_discrete(
        labels = parse(text = grid_label_text[levels(dt_sc$Grid)])
      ) +
      theme_supp(base_size = 14, base_family = "Arial") +
      theme(
        legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          size = 14,
          face = "bold"
        )
      )
    
    ggsave(
      filename = file.path(plot_dir, paste0(file_prefix, "_", sc, ".png")),
      plot = p_sc,
      width = 11,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
}

# ============================================================
# 6. Run R2 and RMSE plots
# ============================================================

# R2 from corr^2
r2_res <- read_metric_files(base_dir, metric = "R2")
plot_metric_by_architecture(r2_res, metric = "R2")

fwrite(
  r2_res,
  file = file.path(plot_dir, "PGS_R2_supplementary_source_data.tsv"),
  sep = "\t"
)

# RMSE
rmse_res <- read_metric_files(base_dir, metric = "RMSE")
plot_metric_by_architecture(rmse_res, metric = "RMSE")

fwrite(
  rmse_res,
  file = file.path(plot_dir, "PGS_RMSE_supplementary_source_data.tsv"),
  sep = "\t"
)

# Combined source data
combined_res <- rbindlist(list(r2_res, rmse_res), fill = TRUE)

fwrite(
  combined_res,
  file = file.path(plot_dir, "PGS_R2_RMSE_supplementary_source_data.tsv"),
  sep = "\t"
)

message("Done. R2 and RMSE supplementary figures saved to: ", plot_dir)

