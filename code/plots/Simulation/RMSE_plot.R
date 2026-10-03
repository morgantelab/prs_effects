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
KEEP_GRID <- "p10e-8_r0.1"

plot_dir <- file.path(OUT_DIR, "plots_RMSE_selected")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Poster-style font sizes
AXIS_TITLE_SIZE <- 34
AXIS_TEXT_SIZE  <- 30
LEGEND_TITLE_SIZE <- 32
LEGEND_TEXT_SIZE  <- 30
PLOT_TITLE_SIZE <- 36

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

theme_rmse <- function(base_family = "Arial") {
  theme_classic(base_family = base_family) +
    theme(
      axis.text.x = element_text(
        size = AXIS_TEXT_SIZE,
        colour = "black",
        angle = 0,
        hjust = 0.5,
        family = base_family
      ),
      axis.text.y = element_text(
        size = AXIS_TEXT_SIZE,
        colour = "black",
        family = base_family
      ),
      axis.title.x = element_text(
        size = AXIS_TITLE_SIZE,
        colour = "black",
        family = base_family
      ),
      axis.title.y = element_text(
        size = AXIS_TITLE_SIZE,
        colour = "black",
        family = base_family
      ),
      legend.title = element_text(
        size = LEGEND_TITLE_SIZE,
        colour = "black",
        family = base_family
      ),
      legend.text = element_text(
        size = LEGEND_TEXT_SIZE,
        colour = "black",
        family = base_family
      ),
      plot.title = element_text(
        size = PLOT_TITLE_SIZE,
        face = "bold",
        hjust = 0.5,
        colour = "black",
        family = base_family
      ),
      axis.line = element_line(linewidth = 1.2, colour = "black"),
      axis.ticks = element_line(linewidth = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )
}

# ============================================================
# 3. Build RMSE results from raw files
# ============================================================
scenario_dirs <- list.dirs(OUT_DIR, recursive = FALSE, full.names = FALSE)
scenario_dirs <- intersect(scenario_dirs, KEEP_SCEN)

res_list <- list()
idx <- 1L

for (scen in scenario_dirs) {
  scen_dir <- file.path(OUT_DIR, scen)
  
  reps <- list.dirs(scen_dir, recursive = FALSE, full.names = FALSE)
  reps <- reps[grepl("^rep", reps)]
  if (length(reps) == 0) next
  
  for (rep in reps) {
    rep_dir <- file.path(scen_dir, rep)
    
    pheno_file <- file.path(rep_dir, "pheno", "simulated_pheno.txt")
    pheno_dt <- read_pheno(pheno_file)
    if (is.null(pheno_dt)) next
    
    gd <- file.path(rep_dir, "grid", KEEP_GRID)
    if (!dir.exists(gd)) next
    
    beta_file <- file.path(gd, "score", "testing_prs_beta_clumped.txt")
    ad_file   <- file.path(gd, "score", "testing_prs_a_d_method2.txt")
    
    prs_beta <- read_prs(beta_file)
    if (!is.null(prs_beta)) {
      dt_beta <- merge(pheno_dt, prs_beta, by = "IID", all = FALSE)
      if (nrow(dt_beta) > 0) {
        res_list[[idx]] <- data.table(
          Scenario = scen,
          Rep      = rep,
          Grid     = KEEP_GRID,
          Type     = "PGS_beta",
          RMSE     = calc_rmse(dt_beta$y_true, dt_beta$y_pred)
        )
        idx <- idx + 1L
      }
    }
    
    prs_ad <- read_prs(ad_file)
    if (!is.null(prs_ad)) {
      dt_ad <- merge(pheno_dt, prs_ad, by = "IID", all = FALSE)
      if (nrow(dt_ad) > 0) {
        res_list[[idx]] <- data.table(
          Scenario = scen,
          Rep      = rep,
          Grid     = KEEP_GRID,
          Type     = "PGS_a_d",
          RMSE     = calc_rmse(dt_ad$y_true, dt_ad$y_pred)
        )
        idx <- idx + 1L
      }
    }
  }
}

if (length(res_list) == 0) {
  stop("No RMSE results were generated. Check file paths and folder structure.")
}

rmse_res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

fwrite(
  rmse_res,
  file.path(plot_dir, "RMSE_results_selected.tsv"),
  sep = "\t"
)

# ============================================================
# 4. Clean labels and factor order
# ============================================================
scen_labels <- c(
  add1000    = "Pure additive",
  mix500_500 = "Combined additive\nand dominance",
  comdom1000 = "Complete dominance"
)

rmse_res <- rmse_res[Scenario %in% KEEP_SCEN & Grid == KEEP_GRID & !is.na(RMSE)]
if (nrow(rmse_res) == 0) {
  stop("After filtering, no RMSE data left.")
}

rmse_res[, Scenario_label := factor(
  Scenario,
  levels = names(scen_labels),
  labels = scen_labels
)]
rmse_res[, Type := factor(Type, levels = c("PGS_a_d", "PGS_beta"))]

# ============================================================
# 5. Colors and labels
# ============================================================
fill_cols <- c(
  "PGS_a_d"  = "#FA7F6F",
  "PGS_beta" = "#82B0D2"
)

type_labels <- c(
  PGS_a_d  = expression(plain(PGS)[italic(a) * "," * italic(d)]),
  PGS_beta = expression(plain(PGS)[bolditalic(beta)])
)

# ============================================================
# 6. Overview plot
# ============================================================
p_rmse <- ggplot(rmse_res, aes(x = Scenario_label, y = RMSE, fill = Type)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.62,
    outlier.shape = 16,
    outlier.size = 2.0,
    outlier.alpha = 0.6,
    linewidth = 0.8
  ) +
  labs(
    x = "Genetic architecture",
    y = "RMSE",
    fill = "PGS type"
  ) +
  scale_fill_manual(
    values = fill_cols,
    labels = type_labels
  ) +
  scale_x_discrete(
    labels = function(x) str_wrap(x, width = 18)
  ) +
  theme_rmse(base_family = "Arial") +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(plot_dir, "PGS_RMSE_overview_p10e-8_r0.1.png"),
  plot = p_rmse,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

message("Done. RMSE plot saved to: ", plot_dir)
