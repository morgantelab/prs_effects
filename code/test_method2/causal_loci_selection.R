#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: causal_loci_selection.R FILTER_DIR out_causal out_ids out_eff SEED N_ADD N_COM N_INCOM")
}

FILTER_DIR <- args[1]
out_causal <- args[2]   # full SNP table
out_ids    <- args[3]   # SNP IDs only (1 col, no header)
out_eff    <- args[4]   # effects table
SEED       <- as.integer(args[5])
N_ADD      <- as.integer(args[6])
N_COM      <- as.integer(args[7])
N_INCOM    <- as.integer(args[8])

# ----- load all available SNPs from filtered BIMs -----
bims <- list.files(FILTER_DIR, pattern = "_filtered\\.bim$", full.names = TRUE)
stopifnot(length(bims) > 0)

all <- rbindlist(lapply(
  bims,
  function(f) fread(f, header = FALSE,
                    col.names = c("chr","snp","cm","pos","a1","a2"))
), use.names = TRUE)

num_loci <- N_ADD + N_COM + N_INCOM
if (num_loci <= 0) stop("Total number of loci (N_ADD + N_COM + N_INCOM) must be > 0.")
if (nrow(all) < num_loci) stop("Not enough SNPs available to sample: needed=", num_loci, " available=", nrow(all))

set.seed(SEED)
sel <- all[sample(.N, num_loci)]

# ----- write causal SNP table & ID list -----
# full table with columns: chr snp cm pos a1 a2
fwrite(sel, out_causal, sep = "\t")

# **FIX**: wrap the vector for fwrite
# one column, no header
fwrite(data.table(sel$snp), out_ids, sep = "\t", col.names = FALSE)

# ----- simulate effects (a, d) with complete / incomplete / additive dominance -----
set.seed(SEED)
a <- rnorm(num_loci, mean = 0, sd = sqrt(2))

# compose and shuffle dominance types
dom_type <- c(
  rep("additive",   N_ADD),
  rep("complete",   N_COM),
  rep("incomplete", N_INCOM)
)
dom_type <- sample(dom_type, num_loci, replace = FALSE)

d <- numeric(num_loci)

# complete dominance: d = ±a
idx_c <- which(dom_type == "complete")
if (length(idx_c) > 0) {
  d[idx_c] <- a[idx_c] * sample(c(-1, 1), length(idx_c), replace = TRUE)
}

# incomplete dominance: d = (±) u * a, u ~ Uniform(0.2, 0.8)
idx_i <- which(dom_type == "incomplete")
if (length(idx_i) > 0) {
  d[idx_i] <- a[idx_i] *
    sample(c(-1, 1), length(idx_i), replace = TRUE) *
    runif(length(idx_i), 0.2, 0.8)
}

# additive remain 0
# idx_a <- which(dom_type == "additive")

true_eff <- data.table(
  rsid = sel$snp,
  a    = a,
  d    = d,
  type = dom_type
)
fwrite(true_eff, out_eff, sep = "\t")
