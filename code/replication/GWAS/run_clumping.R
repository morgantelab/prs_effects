#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
geno_pattern <- args[1]
gwas_pattern <- args[2]
train_ids <- args[3]
out_prefix <- args[4]
plink2_path <- args[5]
p_threshold <- as.numeric(args[6])
ld_thresholds <- args[7]
rep = <- args[8]
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds, ",")))
cat("LD thresholds", ld_thresholds, "\n")

for (r2 in ld_thresholds) {
    r2_str <- format(r2, nsmall = 1)
    cat("Processing replicate", rep, "with LD threshold r² =", r2_str, "\n")

# Load required packages
library(data.table)
chromosomes <- 1:22

for (chr in chromosomes) {
    # Construct file paths for this chromosome
    bfile <- paste0(geno_pattern, chr, "_filtered")
    gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
    for (r2 in ld_thresholds) {
        r2_str <- format(r2, nsmall = 1)
        cat("Processing replicate", rep, "chromosome", chr, "with LD threshold r² =", r2_str, "\n")
        
        # Construct output prefix for this chromosome and threshold
        chr_out_prefix <- paste0(out_prefix, chr, "_r2", r2_str)
        
        # Run clumping
        plink_cmd <- sprintf(
            "%s --bfile %s --clump %s --clump-p1 %s --clump-r2 %s --clump-kb 1000 --keep %s --out %s",
            plink2_path, bfile, gwas_file, p_threshold, r2, train_ids, chr_out_prefix
        )
        system(plink_cmd)
    }
}
