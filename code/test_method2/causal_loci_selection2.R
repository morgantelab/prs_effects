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
out_eff    <- args[4]   # effects table (这里只写 a，d 先置 0，type 可留作占位)
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
fwrite(sel, out_causal, sep = "\t")
fwrite(data.table(sel$snp), out_ids, sep = "\t", col.names = FALSE)

# ----- simulate additive effects only (a); d 留待下一步学习 -----
set.seed(SEED)
a <- rnorm(num_loci, mean = 0, sd = sqrt(2))

true_eff <- data.table(
  rsid = sel$snp,
  a    = a,
  d    = 0,                
  type = "to_learn"        
)
fwrite(true_eff, out_eff, sep = "\t")

cat("Selected ", num_loci, " loci and wrote a (d=0 for now):\n  - ", out_causal, "\n  - ", out_ids, "\n  - ", out_eff, "\n", sep = "")
