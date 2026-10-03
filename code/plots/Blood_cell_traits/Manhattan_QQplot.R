#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

KVIK_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/kvik"

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/Plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c(
  "WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW",
  "Platelet_count", "Plateletcrit", "PDW", "Lymphocyte_perc",
  "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc",
  "Basophill_perc", "Reticulocyte_perc", "MSCV", "HLR_perc"
)

SIG_P <- 5e-8

# =========================
# QQ
# =========================
plot_qq <- function(pvals, out_png, trait) {
  pvals <- pvals[!is.na(pvals) & pvals > 0]
  
  if (length(pvals) == 0) return(NULL)
  
  n <- length(pvals)
  dt <- data.table(
    exp = -log10(ppoints(n)),
    obs = -log10(sort(pvals))
  )
  
  p <- ggplot(dt, aes(exp, obs)) +
    geom_point(color = "#4C72B0", size = 1.2, alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype="dashed", color="red") +
    labs(
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    theme_classic()
  
  ggsave(out_png, p, width = 5, height = 5, dpi = 150)
  cat("✔ QQ saved →", out_png, "\n")
}

# =========================
# Manhattan
# =========================
plot_manhattan <- function(gwas, out_png, trait, sig = 5e-8) {
  
  dt <- gwas[Chromosome %in% 1:22 & !is.na(Wald_P) & !is.na(Basepair)]
  
  chr_len <- dt[, .(max_bp = max(Basepair)), by = Chromosome][order(Chromosome)]
  chr_len[, cum := cumsum(shift(max_bp, fill = 0))]
  
  dt <- merge(dt, chr_len[, .(Chromosome, cum)], by = "Chromosome")
  dt[, pos := Basepair + cum]
  
  ticks <- dt[, .(mid = mean(pos)), by = Chromosome][order(Chromosome)]
  
  p <- ggplot(dt, aes(pos, -log10(Wald_P))) +
    geom_point(color = "#B0B0B0", size = 1, alpha = 0.7) +
    geom_hline(yintercept = -log10(sig), linetype="dashed", color="red") +
    scale_x_continuous(breaks = ticks$mid, labels = ticks$Chromosome) +
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic() +
    theme(panel.grid = element_blank())
  
  ggsave(out_png, p, width = 12, height = 5, dpi = 150)
  cat("✔ Manhattan saved →", out_png, "\n")
}

# =========================
# MAIN
# =========================
for (i in seq_along(TRAITS)) {
  
  trait <- TRAITS[i]
  
  cat("\n=============================\n")
  cat("Processing phenotype", i, "→", trait, "\n")
  
  assoc_file <- file.path(KVIK_DIR, sprintf("kvik_output.step2.pheno%d.assoc", i))
  
  if (!file.exists(assoc_file)) {
    cat("  ❌ Missing:", assoc_file, "\n")
    next
  }
  
  gwas <- fread(assoc_file)
  
  gwas[, Chromosome := as.numeric(Chromosome)]
  gwas[, Basepair := as.numeric(Basepair)]
  gwas[, Wald_P := as.numeric(Wald_P)]
  
  plot_manhattan(
    gwas,
    file.path(OUT_DIR, paste0("manhattan_", trait, ".png")),
    trait
  )
  
  plot_qq(
    gwas$Wald_P,
    file.path(OUT_DIR, paste0("qqplot_", trait, ".png")),
    trait
  )
}

cat("\n🎉 Done. All plots saved to:\n", OUT_DIR, "\n")