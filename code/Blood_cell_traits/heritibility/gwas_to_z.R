#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gwas_to_z.R <in_gwas> <out_z> <mode:additive|dominance>\n",
       "Example: Rscript gwas_to_z.R RDW_ad.txt RDW.add.z.txt additive")
}

in_gwas <- args[1]
out_z   <- args[2]
mode    <- tolower(args[3])

if (!(mode %in% c("additive", "dominance"))) {
  stop("mode must be 'additive' or 'dominance'")
}

dt <- fread(in_gwas, sep = "\t", header = TRUE, data.table = FALSE)

# basic column checks
need_common <- c("SNP", "REF", "ALT")
miss_common <- setdiff(need_common, colnames(dt))
if (length(miss_common) > 0) stop(paste("Missing columns:", paste(miss_common, collapse=", ")))

if (mode == "additive") {
  need <- c("a_hat", "p_a")
  miss <- setdiff(need, colnames(dt))
  if (length(miss) > 0) stop(paste("Missing additive columns:", paste(miss, collapse=", ")))
  
  beta <- dt$a_hat
  pval <- dt$p_a
} else {
  need <- c("d_hat", "p_d")
  miss <- setdiff(need, colnames(dt))
  if (length(miss) > 0) stop(paste("Missing dominance columns:", paste(miss, collapse=", ")))
  
  beta <- dt$d_hat
  pval <- dt$p_d
}

# guard: p cannot be <= 0 (avoid Inf)
pval[pval <= 0 | is.na(pval)] <- .Machine$double.xmin

z <- sign(beta) * sqrt(-2 * log(pval))

out <- data.frame(
  SNP = dt$SNP,
  A1  = dt$ALT,
  A2  = dt$REF,
  Z   = z,
  stringsAsFactors = FALSE
)

dir.create(dirname(out_z), recursive = TRUE, showWarnings = FALSE)
fwrite(out, file = out_z, sep = "\t", quote = FALSE, row.names = FALSE)
