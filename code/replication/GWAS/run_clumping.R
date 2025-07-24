#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
bfile <- args[1]
gwas_file <- args[2]
train_ids <- args[3]
out_prefix <- args[4]
plink2_path <- args[5]
p_threshold <- as.numeric(args[6])
ld_thresh <- args[7]

# Run clumping
cmd <- sprintf(
  "%s --bfile %s --clump %s --clump-p1 %s --clump-r2 %s --clump-kb 1000 --keep %s --out %s",
  plink2_path, bfile, gwas_file, p_threshold, ld_thresh, train_ids, out_prefix
)
system(cmd)