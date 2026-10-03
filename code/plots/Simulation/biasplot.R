suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_centered"

# ---- find PRS prediction files ----
beta_files <- list.files(
  OUT_DIR, pattern = "^testing_prs_beta_clumped\\.txt$",
  recursive = TRUE, full.names = TRUE
)
ad_files <- list.files(
  OUT_DIR, pattern = "^testing_prs_a_d_method2\\.txt$",
  recursive = TRUE, full.names = TRUE
)

all_pred_files <- data.table(
  file = c(beta_files, ad_files),
  Type = c(rep("PRS_beta", length(beta_files)),
           rep("PRS_a_d", length(ad_files)))
)

# ---- parse Scenario / Rep / Grid from path ----
pat <- paste0("^", gsub("([.\\[\\]\\(\\)\\^\\$\\|\\?\\*\\+\\{\\}])", "\\\\\\1", OUT_DIR),
              "/([^/]+)/rep(\\d+)/grid/(p[^/]+)/score/")
m <- str_match(all_pred_files$file, pat)

all_pred_files[, `:=`(
  Scenario = m[,2],
  Rep      = m[,3],
  Grid     = m[,4]
)]

# Split Grid into P and R2
all_pred_files[, P  := str_match(Grid, "p([^_]+)")[,2]]
all_pred_files[, R2 := str_match(Grid, "r(.+)")[,2]]
all_pred_files[, P  := paste0("10E-", sub(".*e-","",P))]
all_pred_files[, xlab2 := paste(P, R2, sep = " • ")]

# ---- helpers ----
read_pheno <- function(scen, rep) {
  f <- file.path(OUT_DIR, scen, paste0("rep", rep), "pheno", "simulated_pheno.txt")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  cn <- tolower(names(dt))
  # common possibilities: FID IID pheno  OR  IID pheno  OR  eid pheno
  fid_col  <- which(cn %in% c("fid","familyid"))
  iid_col  <- which(cn %in% c("iid","eid","id"))
  y_col    <- which(cn %in% c("pheno","phenotype","y_true","y","trait"))
  
  # fallback to positional if needed
  if (length(iid_col) == 0 && ncol(dt) >= 1) iid_col <- 1
  if (length(y_col)   == 0 && ncol(dt) >= 2) y_col   <- 2
  
  res <- data.table(
    IID    = as.character(dt[[iid_col[1]]]),
    y_true = suppressWarnings(as.numeric(dt[[y_col[1]]]))
  )
  if (length(fid_col) > 0) {
    res[, FID := as.character(dt[[fid_col[1]]])]
    setcolorder(res, c("FID","IID","y_true"))
  }
  res[!is.na(y_true)]
}

read_pred <- function(fpath) {
  if (!file.exists(fpath)) return(NULL)
  dt <- tryCatch(fread(fpath), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  cn <- tolower(names(dt))
  # Your case: two columns: IID, PRS
  iid_col <- which(cn %in% c("iid","eid","id"))
  prs_col <- which(cn %in% c("prs","score","y_pred","pred","beta","ad"))
  # fallbacks (positional)
  if (length(iid_col) == 0 && ncol(dt) >= 1) iid_col <- 1
  if (length(prs_col) == 0 && ncol(dt) >= 2) prs_col <- 2
  
  res <- data.table(
    IID    = as.character(dt[[iid_col[1]]]),
    y_pred = suppressWarnings(as.numeric(dt[[prs_col[1]]]))
  )
  # If a FID column exists, keep it
  fid_col <- which(cn %in% c("fid","familyid"))
  if (length(fid_col) > 0) {
    res[, FID := as.character(dt[[fid_col[1]]])]
    setcolorder(res, c("FID","IID","y_pred"))
  }
  res[!is.na(y_pred)]
}

# ---- fit using best available keys (prefer FID+IID, else IID) ----
fit_bias <- function(pheno, pred) {
  if (is.null(pheno) || is.null(pred)) return(list(intercept = NA_real_, slope = NA_real_, n = 0))
  
  common_keys <- intersect(names(pheno), names(pred))
  keys <- intersect(common_keys, c("FID","IID"))
  if (length(keys) >= 2) {
    by_keys <- c("FID","IID")
  } else if ("IID" %in% keys) {
    by_keys <- "IID"
  } else {
    return(list(intercept = NA_real_, slope = NA_real_, n = 0))
  }
  
  dt <- merge(pheno, pred, by = by_keys, allow.cartesian = FALSE)
  dt <- dt[!is.na(y_true) & !is.na(y_pred)]
  if (nrow(dt) == 0) return(list(intercept = NA_real_, slope = NA_real_, n = 0))
  
  fit <- lm(y_true ~ y_pred, data = dt)
  cf  <- coef(fit)
  list(intercept = unname(cf[1]), slope = unname(cf[2]), n = nrow(dt))
}

# ---- compute intercept & slope ----
bias_res <- rbindlist(lapply(seq_len(nrow(all_pred_files)), function(i) {
  row <- all_pred_files[i]
  pheno <- read_pheno(row$Scenario, row$Rep)
  pred  <- read_pred(row$file)
  if (is.null(pheno) || is.null(pred)) {
    return(data.table(
      Scenario=row$Scenario, Rep=row$Rep, P=row$P, R2=row$R2, Type=row$Type,
      Intercept=NA_real_, Slope=NA_real_
    ))
  }
  b <- fit_bias(pheno, pred)
  data.table(
    Scenario=row$Scenario, Rep=row$Rep, P=row$P, R2=row$R2, Type=row$Type,
    Intercept=b$intercept, Slope=b$slope
  )
}), fill=TRUE)

# ---- tidy labels ----
bias_res[, P := factor(P, levels=c("10E-5","10E-8"))]
bias_res[, R2_num := suppressWarnings(as.numeric(R2))]
bias_res[, R2 := factor(R2, levels=unique(R2[order(R2_num)]))]
bias_res[, R2_num := NULL]
bias_res[, xlab2 := paste(P,R2,sep=" • ")]

# save as tab-delimited text
fwrite(bias_res,
       file = file.path(OUT_DIR, "bias_res.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)
# ---- plotting ----
plot_dir <- file.path(OUT_DIR, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
bias_res_path <- file.path(OUT_DIR, "bias_res.txt")
bias_res <- fread(bias_res_path)
# Intercept
p_intercept <- ggplot(bias_res[!is.na(Intercept)],
                      aes(x=xlab2, y=Intercept, fill=Type)) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape=16, outlier.alpha=0.5) +
  facet_wrap(~Scenario, scales="free_x") +
  scale_fill_manual(values=c("PRS_beta"="#377eb8","PRS_a_d"="#e41a1c")) +
  labs(title="Regression intercept (true ~ predicted)",
       x="P-value threshold • R²", y="Intercept", fill="PRS type") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
p_intercept

ggsave(file.path(plot_dir,"PRS_bias_intercept.png"), p_intercept,
       width=14, height=8, dpi=300)

p_intercept <- ggplot(bias_res[!is.na(Intercept)],
                      aes(x = xlab2, y = Intercept, fill = Type)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA   # <-- removes outliers
  ) +
  facet_wrap(~ Scenario, scales = "free_y") +  # <-- facet-specific y-axis
  scale_fill_manual(values = c("PRS_beta" = "#377eb8", "PRS_a_d" = "#e41a1c")) +
  labs(
    title = "Regression intercept: true ~ predicted",
    x = "P-value threshold • R²",
    y = "Intercept",
    fill = "PRS type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )
p_intercept
ggsave(file.path(plot_dir, "PRS_bias_intercept_faceted_freeY.png"),
       p_intercept, width = 14, height = 8, dpi = 300)

# Slope
p_slope <- ggplot(bias_res[!is.na(Slope)],
                  aes(x=xlab2, y=Slope, fill=Type)) +
  geom_boxplot(position=position_dodge(0.8), outlier.shape=16, outlier.alpha=0.5) +
  facet_wrap(~Scenario, scales="free_x") +
  scale_fill_manual(values=c("PRS_beta"="#377eb8","PRS_a_d"="#e41a1c")) +
  labs(title="Regression slope (true ~ predicted)",
       x="P-value threshold • R²", y="Slope", fill="PRS type") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(plot_dir,"PRS_bias_slope.png"), p_slope,
       width=14, height=8, dpi=300)


# --- filter for plotting only ---
bias_res_plot <- copy(bias_res)[R2 == 0.1 & P %in% c(1e-04, 1e-07)]

# label P as "10E-5" and "10E-8"
bias_res_plot[, P_lab := fifelse(abs(P - 1e-04) < .Machine$double.eps^0.5, "10E-5",
                                 fifelse(abs(P - 1e-07) < .Machine$double.eps^0.5, "10E-8", as.character(P)))]
bias_res_plot[, P_lab := factor(P_lab, levels = c("10E-5","10E-8"))]

# ensure red (PRS_a_d) is left and legend matches
bias_res_plot[, Type := factor(Type, levels = c("PRS_a_d","PRS_beta"))]

# use P on x since R2 is fixed at 0.1
# -------- Intercept --------
p_intercept <- ggplot(bias_res_plot[!is.na(Intercept)],
                      aes(x = P_lab, y = Intercept, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_fill_manual(values = c("PRS_a_d" = "#e41a1c", "PRS_beta" = "#377eb8"),
                    name = "PRS type") +
  labs(title = "Regression intercept: true ~ predicted (P ∈ {10E-5, 10E-8}, R² = 0.1)",
       x = "P-value threshold", y = "Intercept") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5), plot.title = element_text(face = "bold"))

p_intercept

ggsave(file.path(plot_dir, "PRS_bias_intercept_P10E-5_10E-8_R0.1.png"),
       p_intercept, width = 12, height = 7, dpi = 300)

# -------- Slope --------
p_slope <- ggplot(bias_res_plot[!is.na(Slope)],
                  aes(x = P_lab, y = Slope, fill = Type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  facet_wrap(~ Scenario, scales = "free_y") +
  scale_fill_manual(values = c("PRS_a_d" = "#e41a1c", "PRS_beta" = "#377eb8"),
                    name = "PRS type") +
  labs(title = "Regression slope: true ~ predicted (P ∈ {10E-5, 10E-8}, R² = 0.1)",
       x = "P-value threshold", y = "Slope") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5), plot.title = element_text(face = "bold"))


p_slope

ggsave(file.path(plot_dir, "PRS_bias_slope_P10E-5_10E-8_R0.1.png"),
       p_slope, width = 12, height = 7, dpi = 300)
