suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---------- Config ----------
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_centered"

# helper functions
read_pheno <- function(path) {
  dt <- fread(path, col.names = c("FID","IID","y_true"))
  dt[, IID := as.character(IID)]
  dt
}

read_prs <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  setnames(dt, c("IID","y_pred"))
  dt[, IID := as.character(IID)]
  dt
}

# ---------- Build MAE + RMSE ----------
scenarios <- list.dirs(OUT_DIR, recursive = FALSE, full.names = FALSE)
res_all <- list()

for (scen in scenarios) {
  reps <- list.dirs(file.path(OUT_DIR, scen), recursive = FALSE, full.names = FALSE)
  reps <- reps[grepl("^rep", reps)]
  
  for (rep in reps) {
    pheno_file <- file.path(OUT_DIR, scen, rep, "pheno", "simulated_pheno.txt")
    if (!file.exists(pheno_file)) next
    pheno <- read_pheno(pheno_file)
    
    grid_dir <- file.path(OUT_DIR, scen, rep, "grid")
    grids <- list.dirs(grid_dir, recursive = FALSE, full.names = TRUE)
    
    for (gd in grids) {
      grid_name <- basename(gd)
      pval <- str_match(grid_name, "p([^_]+)")[,2]
      r2   <- str_match(grid_name, "r(.+)")[,2]
      P <- paste0("10E-", sub(".*e-","",pval))
      R2 <- r2
      
      f_beta <- file.path(gd, "score", "testing_prs_beta_clumped.txt")
      f_ad   <- file.path(gd, "score", "testing_prs_a_d_method2.txt")
      
      # PRS_beta
      if (file.exists(f_beta)) {
        dt <- merge(pheno, read_prs(f_beta), by = "IID", all = FALSE)
        if (nrow(dt) > 0) {
          mae  <- mean(abs(dt$y_true - dt$y_pred))
          rmse <- sqrt(mean((dt$y_true - dt$y_pred)^2))
          res_all[[length(res_all)+1]] <- data.table(
            Scenario=scen, Rep=rep, P=P, R2=R2, Type="PRS_beta",
            MAE=mae, RMSE=rmse
          )
        }
      }
      
      # PRS_a_d
      if (file.exists(f_ad)) {
        dt <- merge(pheno, read_prs(f_ad), by = "IID", all = FALSE)
        if (nrow(dt) > 0) {
          mae  <- mean(abs(dt$y_true - dt$y_pred))
          rmse <- sqrt(mean((dt$y_true - dt$y_pred)^2))
          res_all[[length(res_all)+1]] <- data.table(
            Scenario=scen, Rep=rep, P=P, R2=R2, Type="PRS_a_d",
            MAE=mae, RMSE=rmse
          )
        }
      }
    }
  }
}

acc_res <- rbindlist(res_all, use.names = TRUE, fill = TRUE)
fwrite(acc_res, file.path(OUT_DIR, "MAE_RMSE_results.txt"), sep="\t")

# ---------- Plot MAE ----------
acc_res_path <- file.path(OUT_DIR, "MAE_RMSE_results.txt")
acc_res <- fread(acc_res_path)
acc_res[, P_lab := fifelse(abs(P - 1e-04) < .Machine$double.eps^0.5, "10E-5",
                           fifelse(abs(P - 1e-07) < .Machine$double.eps^0.5, "10E-8", as.character(P)))]
acc_res[, P_lab := factor(P_lab, levels = c("10E-5","10E-8"))]

# Build combined label for x-axis
acc_res[, xlab2 := paste(P_lab, R2, sep = " • ")]

# Reorder PRS types (red = left)
acc_res[, Type := factor(Type, levels = c("PRS_a_d","PRS_beta"))]

plot_dir <- file.path(OUT_DIR, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

p_mae <- ggplot(acc_res, aes(x = xlab2, y = MAE, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  facet_wrap(~Scenario, scales = "free_y") +
  scale_fill_manual(values = c("PRS_a_d" = "#e41a1c", "PRS_beta" = "#377eb8")) +
  labs(title = "Mean Absolute Error (MAE) across scenarios",
       x = "P-value threshold • R²", y = "MAE", fill = "PRS type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

p_mae
ggsave(file.path(plot_dir, "PRS_MAE_boxplot.png"), p_mae, width=14, height=8, dpi=300)

# ---------- Plot RMSE ----------
p_rmse <- ggplot(acc_res, aes(x=xlab2, y=RMSE, fill=Type)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  facet_wrap(~Scenario, scales="free_y") +
  scale_fill_manual(values=c("PRS_beta"="#377eb8","PRS_a_d"="#e41a1c")) +
  labs(title="Root Mean Squared Error (RMSE) across scenarios",
       x="P-value threshold • R²", y="RMSE", fill="PRS type") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.title=element_text(face="bold"))
p_rmse
ggsave(file.path(plot_dir, "PRS_RMSE_boxplot.png"), p_rmse, width=14, height=8, dpi=300)

#-----------Selected P and R2 ---------------------
acc_res <- copy(acc_res)[R2 == 0.1 & P %in% c(1e-04, 1e-07)]
p_mae <- ggplot(acc_res, aes(x=xlab, y=MAE, fill=Type)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  facet_wrap(~Scenario, scales="free_y") +
  scale_fill_manual(values=c("PRS_beta"="#377eb8","PRS_a_d"="#e41a1c")) +
  labs(title="Mean Absolute Error (MAE) across scenarios",
       x="P-value threshold • R²", y="MAE", fill="PRS type") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.title=element_text(face="bold"))
p_mae
ggsave(file.path(plot_dir, "PRS_MAE_boxplot_P10E-5_10E-8_R0.1.png"), p_mae, width=14, height=8, dpi=300)

# ---------- Plot RMSE ----------
p_rmse <- ggplot(acc_res, aes(x=xlab, y=RMSE, fill=Type)) +
  geom_boxplot(position=position_dodge(width=0.8), outlier.shape=NA) +
  facet_wrap(~Scenario, scales="free_y") +
  scale_fill_manual(values=c("PRS_beta"="#377eb8","PRS_a_d"="#e41a1c")) +
  labs(title="Root Mean Squared Error (RMSE) across scenarios",
       x="P-value threshold • R²", y="RMSE", fill="PRS type") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        plot.title=element_text(face="bold"))
p_rmse
ggsave(file.path(plot_dir, "PRS_RMSE_boxplot_P10E-5_10E-8_R0.1.png"), p_rmse, width=14, height=8, dpi=300)


