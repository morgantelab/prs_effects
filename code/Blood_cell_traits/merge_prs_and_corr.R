#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)

score_file <- args[1]   # e.g. score_beta_clumped.sscore
out_prs    <- args[2]
test_ph    <- args[3]
trait      <- args[4]
out_corr   <- args[5]


stopifnot(file.exists(score_file))

# -----------------------------
# Read PRS score
# -----------------------------
dt0 <- fread(score_file)

prs_col <- "SCORE1_SUM"
prs_dt <- dt0[, .(
  IID = as.character(IID),
  PRS = get(prs_col)
)]

fwrite(prs_dt, out_prs, sep = "\t")

# -----------------------------
# Correlation
# -----------------------------
ph_dt <- fread(test_ph)

if (!(trait %in% names(ph_dt))) {
  stop(sprintf("Trait '%s' not found in phenotype file", trait))
}

ph <- ph_dt[, .(
  IID = as.character(IID),
  pheno = get(trait)
)]

ph[, IID := as.character(IID)]

res  <- merge(prs_dt, ph, by = "IID")
rval <- suppressWarnings(cor(res$pheno, res$PRS, use = "complete.obs"))

write.table(
  data.frame(corr_beta = rval, N = nrow(res)),
  out_corr, row.names = FALSE, quote = FALSE, sep = "\t"
)

cat(sprintf("[TEST clumped β] N=%d  cor=%.6f\n", nrow(res), rval))
