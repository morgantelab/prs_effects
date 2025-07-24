#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
effects_file <- args[2]
geno_pattern <- args[3]
keep_file <- args[4]
plink2_path <- args[3]
chromosomes <- as.integer(unlist(strsplit(args[4], ",")))
ld_thresholds <- as.numeric(unlist(strsplit(args[5], ",")))
clump_p <- as.numeric(args[6])
clump_kb <- as.numeric(args[7])
out_prefix <- args[7]
output_file <-args[8]

# Load required packages
library(data.table)

# Read effect estimates
eff <- fread(effects_file)

# Create temporary directory for processing
temp_dir <- file.path(results_dir, sprintf("temp_rep%d", rep))
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Process each LD threshold
for (r2 in ld_thresholds) {
    r2_str <- format(r2, nsmall = 1)
    cat("Processing replicate", rep, "with LD threshold r² =", r2_str, "\n")
    
    clumped_all <- data.table()
    
    # Process each chromosome
    for (chr in chromosomes) {
        cat("  Processing chromosome", chr, "\n")
        
        # Filter effects for current chromosome
        eff_chr <- eff[chromosome == chr]
        if (nrow(eff_chr) == 0) next
        
        # Prepare p-value file for PLINK
        pval_file <- file.path(temp_dir, sprintf("pvals_chr%d.txt", chr))
        fwrite(eff_chr[, .(SNP, P = p_a)], 
               pval_file, 
               sep = "\t", 
               col.names = TRUE,
               na = "NA",
               quote = FALSE)
        
        # Define PLINK inputs/outputs
        plink_file <- sprintf(geno_pattern, chrom = chr)
        out_prefix <- file.path(temp_dir, sprintf("clump_chr%d_r2_%s", chr, r2_str))

        # Build and execute PLINK command
        plink_cmd <- sprintf(
            "%s --bfile %s --keep %s --clump %s --clump-p1 %g --clump-r2 %g --clump-kb %d --out %s",
            plink2_path, plink_file, keep_file,
            pval_file, clump_p, r2, clump_kb, out_prefix
        )
        system(plink_cmd)
        
        # Process clumping results
        clump_file <- paste0(out_prefix, ".clumps")
        clump_chr <- fread(clump_file)
        clump_chr[, CHR := chr]  # Add chromosome column
        clumped_all <- rbindlist(list(clumped_all, clump_chr), fill = TRUE)
    }
    
    # Save final clumped SNPs for this LD threshold
    fwrite(clumped_all, output_file, sep = "\t")
    cat("  Saved", nrow(clumped_all), "clumped SNPs to:", output_file, "\n")
}

# Cleanup temporary files
unlink(temp_dir, recursive = TRUE)
cat("Completed clumping for replicate", rep, "\n")