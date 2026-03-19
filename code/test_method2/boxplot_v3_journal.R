suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal"

# -------------------------
# CONFIG: keep only these
# -------------------------
KEEP_P <- "10E-8"
KEEP_SCEN <- c("add1000", "mix500_500", "comdom1000")

# helper
read_corr_value <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NA_real_)
  suppressWarnings(as.numeric(dt[[1]][1]))
}

# find both beta and ad
files <- list.files(base_dir, pattern="^corr_(beta|ad)\\.txt$", recursive=TRUE, full.names=TRUE)
if (length(files) == 0) stop("No correlation files found. Check base_dir and folder structure.")

all_res <- rbindlist(lapply(files, function(f){
  parts <- str_split(f, "/", simplify=TRUE)
  
  scen <- parts[ncol(parts)-5]
  rep  <- parts[ncol(parts)-4]
  grid <- parts[ncol(parts)-2]  # e.g. p10e-3_r0.1
  
  pval <- str_match(grid, "p([^_]+)")[,2]  # e.g. "10e-8"
  r2   <- str_match(grid, "r(.+)")[,2]
  
  type <- ifelse(str_detect(f, "corr_beta"), "PRS_beta", "PRS_a_d")
  val  <- read_corr_value(f)
  
  data.table(
    Scenario = scen,
    Rep      = rep,
    P        = paste0("10E-", sub(".*e-","",pval)),  # -> "10E-8"
    R2       = r2,
    Type     = type,
    Corr     = val
  )
}), fill = TRUE)

# -------------------------
# FILTER: Scenario + P
# -------------------------
all_res <- all_res[Scenario %in% KEEP_SCEN & P == KEEP_P & R2 == "0.1" & !is.na(Corr)]
if (nrow(all_res) == 0) stop("After filtering, no data left. Check Scenario names and P parsing.")

# -------------------------
# Scenario pretty titles
# -------------------------
scen_labels <- c(
  add1000    = "Pure additive",
  mix500_500 = "Combined additive and dominance",
  comdom1000 = "Complete dominance"
)

all_res[, Scenario_label := factor(Scenario,
                                   levels = names(scen_labels),
                                   labels = scen_labels)]

# -------------------------
# Rename legend labels: PRS -> PGS (DO THIS BEFORE ANY PLOTS)
# -------------------------
# 强制字符，避免旧 factor levels 残留
all_res[, Type := as.character(Type)]

# 明确重命名（最安全）
all_res[Type == "PRS_beta", Type := "PGS_beta"]
all_res[Type == "PRS_a_d",  Type := "PGS_a_d"]

# 再设为有序因子（控制图例顺序）
all_res[, Type := factor(Type, levels = c("PGS_a_d", "PGS_beta"))]

# consistent blue/orange like your example
fill_cols <- c("PGS_beta" = "#82B0D2",
               "PGS_a_d"  = "#FA7F6F")

# output dir
plot_dir <- file.path(base_dir, "plots_10E-8_selected_scen")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# (A) one PNG per Scenario
# -------------------------
for (sc in KEEP_SCEN) {
  dt <- all_res[Scenario == sc]
  if (nrow(dt) == 0) next
  
  p <- ggplot(dt, aes(x = xlab, y = Corr, fill = Type)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.shape = 16, outlier.alpha = 0.5) +
    labs(
      title = scen_labels[sc],
      x = "R²",
      y = "Correlation",
      fill = "PGS type"
    ) +
    scale_fill_manual(values = fill_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(plot_dir, sprintf("PGS_corr_boxplot_%s.png", sc)),
         p, width = 9, height = 6, dpi = 300)
}

# -------------------------
# (B) one faceted overview
# -------------------------
p_all <- ggplot(all_res, aes(x = Type, y = Corr, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Scenario_label, scales = "free_x", nrow = 1) +
  labs(
    title = "PGS correlation across genetic architectures",
    x = "R²",
    y = "Correlation",
    fill = "PGS type"
  ) +
  scale_fill_manual(values = fill_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(file.path(plot_dir, "PGS_corr_boxplot_faceted_selected.png"),
       p_all, width = 14, height = 5.5, dpi = 300)
