#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# User configuration
# ============================================================
KVIK_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/kvik"

# save all plots here
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/Plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

detect_traits <- function(kvik_dir) {
  d <- list.dirs(kvik_dir, full.names = FALSE, recursive = FALSE)
  d <- d[!is.na(d) & d != ""]
  d <- setdiff(d, c("plots"))
  d
}

TRAITS_FALLBACK <- c("HEIGHT", "BMI")

TRAITS <- detect_traits(KVIK_DIR)
if (length(TRAITS) == 0L) {
  TRAITS <- TRAITS_FALLBACK
}

SIG_P <- 5e-8

# ============================================================
# Helper: QQ plot
# ============================================================
plot_qq <- function(pvals, out_png) {
  
  pvals <- pvals[!is.na(pvals) & pvals > 0]
  
  if (length(pvals) == 0) {
    cat("  ❌ No valid p-values for QQ plot → skip\n")
    return(NULL)
  }
  
  n <- length(pvals)
  dt <- data.table(
    exp = -log10(ppoints(n)),
    obs = -log10(sort(pvals))
  )
  
  p <- ggplot(dt, aes(x = exp, y = obs)) +
    geom_point(color = "#4D79A6", size = 1.2, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  
  ggsave(out_png, p, width = 5, height = 5, dpi = 150)
  cat("✔ QQ saved →", out_png, "\n")
}

# ============================================================
# Helper: Manhattan plot
# ============================================================
plot_manhattan <- function(gwas, out_png, sig = 5e-8) {
  
  dt <- copy(gwas)
  dt <- dt[Chromosome %in% 1:22 & !is.na(Wald_P) & !is.na(Basepair)]
  dt <- dt[Wald_P > 0]
  
  if (nrow(dt) == 0) {
    cat("  ❌ No valid rows for Manhattan plot → skip\n")
    return(NULL)
  }
  
  chr_lengths <- dt[, .(max_bp = max(Basepair, na.rm = TRUE)), by = Chromosome][order(Chromosome)]
  chr_lengths[, cum_start := cumsum(shift(as.numeric(max_bp), fill = 0))]
  
  dt <- merge(dt, chr_lengths[, .(Chromosome, cum_start)], by = "Chromosome")
  dt[, pos_cum := Basepair + cum_start]
  
  chr_span <- dt[, .(
    chr_start = min(pos_cum, na.rm = TRUE),
    chr_end   = max(pos_cum, na.rm = TRUE)
  ), by = Chromosome][order(Chromosome)]
  chr_span[, mid := (chr_start + chr_end) / 2]
  
  sig_line <- -log10(sig)
  
  p <- ggplot(dt, aes(x = pos_cum, y = -log10(Wald_P))) +
    geom_point(color = "#4D79A6", size = 0.8, alpha = 0.75) +
    geom_hline(yintercept = sig_line, linetype = "dashed", color = "red") +
    scale_x_continuous(
      breaks = chr_span$mid,
      labels = chr_span$Chromosome
    ) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid = element_blank()
    )
  
  ggsave(out_png, p, width = 12, height = 5, dpi = 150)
  cat("✔ Manhattan saved →", out_png, "\n")
}

# ============================================================
# Main — iterate traits
# ============================================================
for (trait in TRAITS) {
  
  cat("\n=============================\n")
  cat("Processing trait →", trait, "\n")
  
  assoc_file <- file.path(KVIK_DIR, trait, "kvik_output.step2.assoc")
  
  if (!file.exists(assoc_file)) {
    cat("  ❌ assoc file not found → skip: ", assoc_file, "\n")
    next
  }
  
  gwas <- fread(assoc_file)
  
  required_cols <- c("Chromosome", "Basepair", "Wald_P")
  if (!all(required_cols %in% names(gwas))) {
    stop(
      "Missing columns in assoc for trait: ", trait,
      "\nNeed: ", paste(required_cols, collapse = ", "),
      "\nHave: ", paste(names(gwas), collapse = ", ")
    )
  }
  
  gwas[, Chromosome := as.numeric(Chromosome)]
  gwas[, Basepair := as.numeric(Basepair)]
  gwas[, Wald_P := as.numeric(Wald_P)]
  
  manhattan_png <- file.path(OUT_DIR, paste0("manhattan_", trait, ".png"))
  plot_manhattan(gwas, manhattan_png, sig = SIG_P)
  
  qq_png <- file.path(OUT_DIR, paste0("qqplot_", trait, ".png"))
  plot_qq(gwas$Wald_P, qq_png)
}

cat("\n🎉 All traits processed.\n")
cat("Plots saved to: ", OUT_DIR, "\n")