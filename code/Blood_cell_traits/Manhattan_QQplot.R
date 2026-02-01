#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# User configuration
# ============================================================

KVIK_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/kvik"

TRAITS <- c(
  "WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW",
  "Platelet_count", "Plateletcrit", "PDW", "Lymphocyte_perc",
  "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc",
  "Basophill_perc", "Reticulocyte_perc", "MSCV", "HLR_perc"
)

SIG_P <- 5e-8   # genome-wide threshold


# ============================================================
# Helper: QQ plot
# ============================================================
plot_qq <- function(pvals, out_png, trait) {
  
  pvals <- pvals[!is.na(pvals) & pvals > 0]
  
  n <- length(pvals)
  exp <- -log10(ppoints(n))
  obs <- -log10(sort(pvals))
  
  dt <- data.table(exp = exp, obs = obs)
  
  p <- ggplot(dt, aes(x = exp, y = obs)) +
    geom_point(color = "#4C72B0", size = 1.2, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype="dashed") +
    labs(
      title = paste("QQ Plot —", trait),
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    theme_classic()
  
  ggsave(out_png, p, width = 5, height = 5, dpi = 150)
  cat("✔ QQ saved →", out_png, "\n")
}


# ============================================================
# Helper: Manhattan plot
# ============================================================
plot_manhattan <- function(gwas, out_png, trait, sig = 5e-8) {
  
  dt <- gwas[Chromosome %in% 1:22 & !is.na(Wald_P)]
  
  # cumulative positions
  chr_lengths <- dt[, .(max_bp = max(Basepair)), by = Chromosome][order(Chromosome)]
  chr_lengths[, cum_start := cumsum(shift(max_bp, fill = 0))]
  
  dt <- merge(dt, chr_lengths[, .(Chromosome, cum_start)], by = "Chromosome")
  dt[, pos_cum := Basepair + cum_start]
  
  ticks <- dt[, .(mid = mean(pos_cum)), by = Chromosome][order(Chromosome)]
  
  sig_line <- -log10(sig)
  
  p <- ggplot(dt, aes(x = pos_cum, y = -log10(Wald_P))) +
    geom_point(color = "#B0B0B0", size = 1, alpha = 0.7) +
    geom_hline(yintercept = sig_line, linetype = "dashed", color = "red") +
    scale_x_continuous(
      breaks = ticks$mid,
      labels = ticks$Chromosome
    ) +
    labs(
      title = paste("Manhattan Plot —", trait),
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave(out_png, p, width = 12, height = 5, dpi = 150)
  cat("✔ Manhattan saved →", out_png, "\n")
}


# ============================================================
# Main — iterate through pheno1…pheno16
# ============================================================

for (i in seq_along(TRAITS)) {
  
  trait <- TRAITS[i]
  
  cat("\n=============================\n")
  cat("Processing phenotype", i, "→", trait, "\n")
  
  assoc_file <- file.path(KVIK_DIR, sprintf("kvik_output.step2.pheno%d.assoc", i))
  
  if (!file.exists(assoc_file)) {
    cat("  ❌ assoc file not found → skip\n")
    next
  }
  
  # Output directory
  outdir <- file.path(KVIK_DIR, "plots", trait)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Load GWAS
  gwas <- fread(assoc_file)
  
  # Ensure needed columns exist
  required_cols <- c("Chromosome", "Basepair", "Wald_P")
  if (!all(required_cols %in% names(gwas))) {
    stop("Missing columns in assoc for phenotype ", i)
  }
  
  # Convert types
  gwas[, Chromosome := as.numeric(Chromosome)]
  gwas[, Basepair := as.numeric(Basepair)]
  gwas[, Wald_P := as.numeric(Wald_P)]
  
  # Manhattan
  manhattan_png <- file.path(outdir, paste0("manhattan_", trait, ".png"))
  plot_manhattan(gwas, manhattan_png, trait, sig = SIG_P)
  
  # QQ
  qq_png <- file.path(outdir, paste0("qqplot_", trait, ".png"))
  plot_qq(gwas$Wald_P, qq_png, trait)
}

cat("\n🎉 All phenotypes processed.\n")
