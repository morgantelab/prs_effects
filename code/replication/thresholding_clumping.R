#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
effects_file <- args[2]
geno_pattern <- args[3]
keep_file <- args[4]
plink2_path <- args[5]
ld_thresholds_str <- args[6]
clump_p <- as.numeric(args[7])
clump_kb <- as.numeric(args[8])
out_dir <- args[9]  # Output directory for p-threshold

cat("\nConverted values:\n")
cat("Effect Files:", effects_file, "\n")
cat("rep:", rep, "\n")
cat("Geno File Names:", geno_pattern, "\n")
cat("Train ID file:", keep_file, "\n")
cat("clump p value:", clump_p, "\n")
cat("clump kb:", clump_kb, "\n")
cat("Output directory:", out_dir, "\n")
cat("LD_thresholds:", ld_thresholds_str, "\n")

# Load required packages
library(data.table)

# Read effect estimates
eff <- fread(effects_file)

# Create temporary directory for processing
temp_dir <- file.path(out_dir, sprintf("temp_rep%d", rep))
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Convert LD thresholds to numeric vector
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))
cat("LD thresholds", ld_thresholds, "\n")

# Process each LD threshold
for (r2 in ld_thresholds) {
    r2_str <- format(r2, nsmall = 1)
    cat("Processing replicate", rep, "with LD threshold r² =", r2_str, "\n")
    
    clumped_all <- data.table()
    
    # Process each chromosome
    for (chr in 1:22) {
        cat("  Processing chromosome", chr, "\n")
        
        # Filter effects for current chromosome
        eff_chr <- eff[eff$chromosome == chr, ]
        if (nrow(eff_chr) == 0) next
        
        # Prepare p-value file for PLINK
        pval_file <- file.path(temp_dir, sprintf("pvals_chr%d.txt", chr))
        pval_data <- data.frame(SNP = eff_chr$SNP, P = eff_chr$p_a)
        fwrite(pval_data,  # Use the new data frame
            pval_file, 
            sep = "\t", 
            col.names = TRUE,
            na = "NA",
            quote = FALSE)
        
        # Define PLINK inputs/outputs
        plink_file <- paste0(geno_pattern, chr, "_filtered")
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
        if (file.exists(clump_file)) {
            clump_chr <- fread(clump_file)
            clump_chr[, CHR := chr]  # Add chromosome column
            clumped_all <- rbindlist(list(clumped_all, clump_chr), fill = TRUE)
        }
    }
    
    # Save final clumped SNPs for this LD threshold
    output_file <- file.path(out_dir, sprintf("clumped_snps_rep%d_r2_%.1f.txt", rep, r2))
    fwrite(clumped_all, output_file, sep = "\t")
    cat("  Saved", nrow(clumped_all), "clumped SNPs to:", output_file, "\n")
}

# Cleanup temporary files
unlink(temp_dir, recursive = TRUE)
cat("Completed clumping for replicate", rep, "\n")