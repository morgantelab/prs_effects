suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2"

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
  
  pval <- str_match(grid, "p([^_]+)")[,2]
  r2   <- str_match(grid, "r(.+)")[,2]
  
  type <- ifelse(str_detect(f, "corr_beta"), "PRS_beta", "PRS_a_d")
  val  <- read_corr_value(f)
  
  data.table(
    Scenario=scen,
    Rep=rep,
    P=paste0("10E-", sub(".*e-","", pval)),  # 10E-3, 10E-5, 10E-8
    R2=r2,
    Type=type,
    Corr=val
  )
}))

# ---- FILTER: only Scenario=add1000, P in {10E-5, 10E-8}, R2 = 0.1 ----
dt <- all_res[
  Scenario == "add1000" &
    P %in% c("10E-5", "10E-8") &
    R2 %in% c("0.1", "0.10") &
    !is.na(Corr)
]

if (nrow(dt) == 0) stop("No rows after filtering (add1000, P in {10E-5,10E-8}, R2=0.1).")

# tidy factors for clean ordering
dt[, P := factor(P, levels = c("10E-5", "10E-8"))]
dt[, R2 := "0.1"]  # normalize label

# single plot, no facet
plot_dir <- file.path(base_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(dt, aes(x = P, y = Corr, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = 16, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("PRS_beta" = "#377eb8", "PRS_a_d" = "#e41a1c")) +
  labs(
    title = "PRS correlation with phenotype",
    x = "P-value threshold",
    y = "Correlation",
    fill = "PRS type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(plot_dir, "PRS_corr_add1000_P_10E-5_10E-8_R2_0.1.png"),
       p, width = 8, height = 6, dpi = 300)





