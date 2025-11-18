suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# =======================
# Config
# =======================
OUT_DIR      <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V3"
N_REPS_WANT  <- 5           # choose 3 or 4
SEED         <- 2025
STRATA_MODE  <- "decile"  # or "decile"

# =======================
# Helpers
# =======================

# Read phenotype
read_pheno <- function(scen, rep) {
  f <- file.path(OUT_DIR, scen, paste0("rep", rep), "pheno", "simulated_pheno.txt")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  if (ncol(dt) >= 3) {
    setnames(dt, old = names(dt)[1:3], new = c("FID","IID","y_true"))
  } else {
    stop("Unexpected format in ", f)
  }
  dt[, IID := as.character(IID)]
  dt[, y_true := as.numeric(y_true)]
  dt[!is.na(y_true)]
}

# Read PRS (IID, score)
read_prs_two_cols <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  setnames(dt, old = names(dt)[1:2], new = c("IID","y_pred"))
  dt[, IID := as.character(IID)]
  dt[, y_pred := as.numeric(y_pred)]
  dt[!is.na(y_pred)]
}

# Stratify into quintiles/deciles
assign_strata <- function(y, mode = c("quintile","decile")) {
  mode <- match.arg(mode)
  k <- if (mode == "quintile") 5 else 10
  brks <- quantile(y, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  brks <- unique(brks)
  if (length(brks) <= 2) {
    rk <- frank(y, ties.method = "average", na.last = "keep")
    cut(rk,
        breaks = quantile(rk, probs = seq(0,1,length.out=k+1), na.rm=TRUE),
        include.lowest = TRUE,
        labels = paste0(if (k==5) "Q" else "D", 1:k))
  } else {
    cut(y, breaks = brks, include.lowest = TRUE,
        labels = paste0(if (k==5) "Q" else "D", 1:k))
  }
}

# Parse P and R2
parse_grid_labels <- function(grid_name) {
  p_raw <- str_match(grid_name, "p([^_]+)")[,2]
  r2    <- str_match(grid_name, "r(.+)")[,2]
  list(P = paste0("10E-", sub(".*e-","", p_raw)), R2 = r2)
}

# =======================
# Step 1. Pick same reps
# =======================
set.seed(SEED)
all_reps <- as.character(1:20)
picked_reps <- sort(sample(all_reps, size = N_REPS_WANT))
message("Picked reps: ", paste(paste0("rep", picked_reps), collapse = ", "))

scenarios <- list.dirs(OUT_DIR, recursive = FALSE, full.names = FALSE)
picked <- data.table(
  Scenario = rep(scenarios, each = length(picked_reps)),
  Rep      = rep(picked_reps, times = length(scenarios))
)

# =======================
# Step 2. Build abs error table
# =======================
abs_res <- rbindlist(lapply(1:nrow(picked), function(i) {
  scen <- picked$Scenario[i]; rep <- picked$Rep[i]
  pheno <- read_pheno(scen, rep)
  if (is.null(pheno)) return(NULL)
  
  grid_root <- file.path(OUT_DIR, scen, paste0("rep", rep), "grid")
  if (!dir.exists(grid_root)) return(NULL)
  grids <- list.dirs(grid_root, recursive = FALSE, full.names = TRUE)
  if (length(grids) == 0) return(NULL)
  
  rbindlist(lapply(grids, function(gd) {
    grid_name <- basename(gd)
    labs <- parse_grid_labels(grid_name)
    
    f_beta <- file.path(gd, "score", "testing_prs_beta_clumped.txt")
    f_ad   <- file.path(gd, "score", "testing_prs_a_d_method2.txt")
    
    out <- list()
    
    if (file.exists(f_beta)) {
      prs <- read_prs_two_cols(f_beta)
      if (!is.null(prs)) {
        dt <- merge(pheno, prs, by = "IID")
        if (nrow(dt)) {
          dt[, Stratum := assign_strata(y_true, STRATA_MODE)]
          out[[length(out)+1]] <- dt[, .(Scenario = scen, Rep = rep,
                                         Grid = grid_name, P = labs$P, R2 = labs$R2,
                                         Type = "PRS_beta",
                                         Stratum,
                                         AbsErr = abs(y_true - y_pred))]
        }
      }
    }
    if (file.exists(f_ad)) {
      prs <- read_prs_two_cols(f_ad)
      if (!is.null(prs)) {
        dt <- merge(pheno, prs, by = "IID")
        if (nrow(dt)) {
          dt[, Stratum := assign_strata(y_true, STRATA_MODE)]
          out[[length(out)+1]] <- dt[, .(Scenario = scen, Rep = rep,
                                         Grid = grid_name, P = labs$P, R2 = labs$R2,
                                         Type = "PRS_a_d",
                                         Stratum,
                                         AbsErr = abs(y_true - y_pred))]
        }
      }
    }
    rbindlist(out, use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

stopifnot(nrow(abs_res) > 0)

# =======================
# Step 3. Clean factors
# =======================
if (STRATA_MODE == "quintile") {
  abs_res[, Stratum := factor(Stratum, levels = paste0("Q", 1:5))]
} else {
  abs_res[, Stratum := factor(Stratum, levels = paste0("D", 1:10))]
}
abs_res[, Type := factor(Type, levels = c("PRS_beta","PRS_a_d"))]
rep_levels <- paste0("rep", picked_reps)
abs_res[, Rep := factor(paste0("rep", Rep), levels = rep_levels)]

# =======================
# Step 4. Plot
# =======================
plot_dir <- file.path(OUT_DIR, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

abs_res[, Type := factor(Type, levels = c("PRS_a_d","PRS_beta"))]

p_abs <- ggplot(abs_res, aes(x = Stratum, y = AbsErr, fill = Type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  facet_grid(Scenario ~ Rep, scales = "free_y", drop = FALSE) +
  scale_fill_manual(values = c("PRS_a_d"="#e41a1c","PRS_beta"="#377eb8"),
                    name = "PRS method", drop = FALSE) +
  labs(title = "Absolute error",
       x = paste0(if (STRATA_MODE=="quintile") "Quintile" else "Decile",
                  " (by true phenotype)"),
       y = "Absolute error |true − PRS|") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title  = element_text(face = "bold"))

print(p_abs)

ggsave(file.path(plot_dir,
                 sprintf("abs_error_by_%s_same_reps.png", STRATA_MODE)),
       p_abs, width = 16, height = 9, dpi = 300)

# =======================
# Step 5. Save table
# =======================
fwrite(abs_res,
       file.path(OUT_DIR,
                 sprintf("abs_error_table_%s_same_reps.txt", STRATA_MODE)),
       sep = "\t", quote = FALSE)
