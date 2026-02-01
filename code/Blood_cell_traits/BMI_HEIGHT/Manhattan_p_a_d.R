library(data.table)
library(ggplot2)
library(ggrepel)

# Base OUT_DIR for HEIGHT_BMI
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"

# Two phenotypes, handled separately
PHENOS <- c("HEIGHT", "BMI")

# If you want different trait names from pheno folder name, change here:
PHENO_TO_TRAIT <- list(
  HEIGHT = "HEIGHT",
  BMI    = "BMI"
)

# Helper: locate effect dir for a given pheno
get_effect_dir <- function(out_dir, pheno) {
  d1 <- file.path(out_dir, pheno, "effects")
  d2 <- file.path(out_dir, "effects")
  if (dir.exists(d1)) return(d1)
  if (dir.exists(d2)) return(d2)
  stop("Cannot find effects directory for pheno = ", pheno,
       "\nTried: ", d1, "\n   and: ", d2)
}

manhattan_plot <- function(dt, p_col, trait, out_png) {
  setDT(dt)
  
  # Required columns
  need_cols <- c("CHR", "POS", "SNP", p_col)
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  # ---- sanitize types (avoid integer overflow / NA join issues) ----
  dt[, CHR := as.integer(as.character(CHR))]
  dt[, POS := as.numeric(as.character(POS))]
  dt[, pval := as.numeric(as.character(get(p_col)))]
  
  dt <- dt[CHR %in% 1:22]
  dt <- dt[!is.na(CHR) & !is.na(POS) & !is.na(pval)]
  
  # guard p==0 -> use smallest positive
  dt[pval <= 0, pval := .Machine$double.xmin]
  dt[!is.finite(pval), pval := NA_real_]
  dt <- dt[!is.na(pval)]
  
  # ---- chr lengths + cumulative starts (numeric cumsum) ----
  chr_info <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR][order(CHR)]
  chr_info[, cum_start := cumsum(shift(as.numeric(chr_len), fill = 0))]
  
  # join cum_start back
  dt[chr_info, cum_start := i.cum_start, on = "CHR"]
  dt[, pos_cum := POS + cum_start]
  
  # x-axis midpoints by span (not mean)
  chr_span <- dt[, .(
    chr_start = min(pos_cum, na.rm = TRUE),
    chr_end   = max(pos_cum, na.rm = TRUE)
  ), by = CHR][order(CHR)]
  chr_span[, mid := (chr_start + chr_end) / 2]
  
  # ---- thresholds ----
  n_snps <- nrow(dt)
  bonf_line <- -log10(0.05 / n_snps)
  gw_line   <- -log10(5e-8)
  
  # ---- -log10(p) ----
  dt[, logp := -log10(pval)]
  dt[, sig_bonf := pval < (0.05 / n_snps)]
  
  p <- ggplot(dt, aes(x = pos_cum, y = logp)) +
    geom_point(size = 0.8, alpha = 0.75) +
    geom_hline(yintercept = bonf_line, linetype = "dashed") +
    geom_hline(yintercept = gw_line,   linetype = "dashed") +
    scale_x_continuous(breaks = chr_span$mid, labels = chr_span$CHR) +
    labs(
      title = paste0("Manhattan Plot - ", trait, " (", p_col, ")"),
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.title  = element_text(size = 12)
    )
  
  # label: one SNP per chr among Bonferroni-significant SNPs
  label_dt <- dt[sig_bonf == TRUE, .SD[which.min(pval)], by = CHR]
  
  if (nrow(label_dt) > 0) {
    p <- p +
      geom_text_repel(
        data = label_dt,
        aes(label = SNP),
        size = 3,
        max.overlaps = Inf
      )
  }
  
  ggsave(out_png, p, width = 12, height = 5, dpi = 150)
  cat("Saved:", out_png, "\n")
}

# ======================================================
# MAIN LOOP (per pheno)
# ======================================================
for (pheno in PHENOS) {
  
  effect_dir <- get_effect_dir(OUT_DIR, pheno)
  trait <- PHENO_TO_TRAIT[[pheno]]
  
  # default file name: <trait>_ad.txt
  file_path <- file.path(effect_dir, paste0(trait, "_ad.txt"))
  
  if (!file.exists(file_path)) {
    cat("Missing:", file_path, "\n")
    next
  }
  
  dt <- fread(file_path)
  
  if (!("p_a" %in% names(dt)) || !("p_d" %in% names(dt))) {
    stop("File missing p_a or p_d: ", file_path)
  }
  
  plotdir <- file.path(effect_dir, "manhattan_plots")
  dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pa <- file.path(plotdir, paste0(trait, "_manhattan_p_a.png"))
  out_pd <- file.path(plotdir, paste0(trait, "_manhattan_p_d.png"))
  
  manhattan_plot(dt, "p_a", trait, out_pa)
  manhattan_plot(dt, "p_d", trait, out_pd)
}

cat("Done.\n")
