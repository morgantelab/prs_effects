suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---------- Config ----------
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal"

KEEP_P    <- "10E-8"
KEEP_R2   <- "0.1"
KEEP_SCEN <- c("add1000", "mix500_500", "comdom1000")

# ---------- Load results ----------
acc_res_path <- file.path(OUT_DIR, "MAE_RMSE_results.txt")
if (!file.exists(acc_res_path)) stop("MAE_RMSE_results.txt not found: ", acc_res_path)
acc_res <- fread(acc_res_path)

# ---------- Make sure P / R2 are comparable as strings ----------
acc_res[, P_chr  := as.character(P)]
acc_res[, R2_chr := as.character(R2)]

# If P is numeric (e.g., 1e-07), map to 10E-8 label
suppressWarnings(acc_res[, P_num := as.numeric(P_chr)])
acc_res[, P_lab := P_chr]
acc_res[!is.na(P_num) & abs(P_num - 1e-07) < .Machine$double.eps^0.5, P_lab := "10E-8"]

# ---------- Filter: only P=10E-8, R2=0.1, and the 3 scenarios ----------
acc_res <- acc_res[
  Scenario %in% KEEP_SCEN &
    P_lab == KEEP_P &
    R2_chr == KEEP_R2
]

if (nrow(acc_res) == 0) {
  stop("No rows left after filtering. Check values in MAE_RMSE_results.txt")
}

# ---------- Rename PRS -> PGS labels (exactly as requested) ----------
acc_res[, Type := as.character(Type)]
acc_res[Type == "PRS_beta", Type := "PGS_beta"]
acc_res[Type == "PRS_a_d",  Type := "PGS_a_d"]
acc_res[, Type := factor(Type, levels = c("PGS_a_d", "PGS_beta"))]

# ---------- Colors (blue/orange) ----------
fill_cols <- c(
  "PGS_beta" = "#82B0D2",
  "PGS_a_d"  = "#FA7F6F"
)

# ---------- Plot RMSE ----------
plot_dir <- file.path(base_dir, "plots_10E-8_selected_scen")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Since P and R2 are fixed, x-axis should be PGS type (otherwise only one x level)
p_rmse <- ggplot(acc_res, aes(x = Type, y = RMSE, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, width = 0.65) +
  facet_wrap(~Scenario, scales = "free_y") +
  scale_fill_manual(values = fill_cols) +
  labs(
    title = "Root Mean Squared Error (RMSE) across scenarios",
    x = NULL,
    y = "RMSE",
    fill = "PGS type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_rmse)

ggsave(
  file.path(plot_dir, "PGS_RMSE_boxplot_P10E-8_R0.1_3scenarios.png"),
  p_rmse,
  width = 12,
  height = 6,
  dpi = 300
)