suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V3"

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
  # break path into components
  parts <- str_split(f, "/", simplify=TRUE)
  scen <- parts[ncol(parts)-5]
  rep  <- parts[ncol(parts)-4]
  grid <- parts[ncol(parts)-2]  # e.g. p10e-3_r0.1
  
  # extract P and R2
  pval <- str_match(grid, "p([^_]+)")[,2]
  r2   <- str_match(grid, "r(.+)")[,2]
  
  type <- ifelse(str_detect(f, "corr_beta"), "PRS_beta", "PRS_a_d")
  val  <- read_corr_value(f)
  
  data.table(
    Scenario=scen,
    Rep=rep,
    P=paste0("10E-", sub(".*e-","",pval)),  # keep 10E-3, 10E-5, 10E-8 style
    R2=r2,
    Type=type,
    Corr=val
  )
}))

# create combined label for x-axis
all_res[, xlab := paste(Scenario, P, R2)]

ggplot(all_res, aes(x=xlab, y=Corr, fill=Type)) +
  geom_boxplot(position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("PRS_beta"="#377eb8", "PRS_a_d"="#e41a1c")) +
  labs(title="PRS correlation with phenotype across scenarios, P, R2",
       x="Scenario + P + R2", y="Correlation") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(base_dir, "PRS_corr_boxplot.png"), width=14, height=6)



# ---- tidy labels & ordering (optional but helpful) ----
# order P as 10E-3 < 10E-5 < 10E-8
all_res[, P := factor(P, levels = c("10E-3","10E-5","10E-8"))]

# make R2 a nicely ordered factor (e.g., 0.05, 0.1, 0.2 ...)
# keep original string but order numerically
all_res[, R2_num := suppressWarnings(as.numeric(R2))]
all_res[, R2 := factor(R2, levels = unique(R2[order(R2_num)]))]
all_res[, R2_num := NULL]

# x-axis label within each scenario: just P + R2 (Scenario is fixed per plot)
all_res[, xlab2 := paste(P, R2, sep = " • ")]

# output dir for plots
plot_dir <- file.path(base_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---- (A) one PNG per Scenario ----
unique_scen <- sort(unique(all_res$Scenario))

for (sc in unique_scen) {
  dt <- all_res[Scenario == sc & !is.na(Corr)]
  if (nrow(dt) == 0) next
  
  p <- ggplot(dt, aes(x = xlab2, y = Corr, fill = Type)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = 16, outlier.alpha = 0.5) +
    scale_fill_manual(values = c("PRS_beta" = "#377eb8", "PRS_a_d" = "#e41a1c")) +
    labs(
      title = sprintf("correlation with phenotype", sc),
      x = "P-value threshold • R²",
      y = "Correlation",
      fill = "PRS type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  fout <- file.path(plot_dir, sprintf("PRS_corr_boxplot_%s.png", sc))
  ggsave(fout, p, width = 10, height = 6, dpi = 300)
}

# ---- (B) an overview faceted figure (all scenarios together) ----
p_all <- ggplot(all_res[!is.na(Corr)], aes(x = xlab2, y = Corr, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("PRS_beta" = "#377eb8", "PRS_a_d" = "#e41a1c")) +
  facet_wrap(~ Scenario, scales = "free_x") +
  labs(
    title = "PRS correlation across scenarios (faceted)",
    x = "P-value threshold • R²",
    y = "Correlation",
    fill = "PRS type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(plot_dir, "PRS_corr_boxplot_all_scenarios_faceted.png"),
       p_all, width = 14, height = 8, dpi = 300)

