#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
geno_pattern <- args[1]        # Genotype file pattern
gwas_pattern <- args[2]        # GWAS results pattern
clump_pattern <- args[3]       # Clumping results pattern
test_ids_file <- args[4]       # File with test IDs
output_prefix <- args[5]       # Output prefix
plink2_path <- args[6]         # Path to PLINK2
ld_thresholds_str <- args[7]   # Comma-separated LD thresholds
results_dir <- args[8]         # Results directory
rep <- args[9]                # Replicate number

# Convert LD thresholds to numeric vector
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))
cat("LD thresholds:", ld_thresholds, "\n")

# Create temporary directory
temp_dir <- file.path(results_dir, "temp_prs")
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
  cat("Created temporary directory:", temp_dir, "\n")
}

# Loading required packages
library(data.table)
library(genio)

# Load test IDs
test_ids <- fread(test_ids_file, header = FALSE)[[1]]
cat("Loaded", length(test_ids), "test IDs\n")

# Initialize PRS matrix
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = length(ld_thresholds))
colnames(prs_matrix) <- sprintf("r2_%s", format(ld_thresholds, nsmall = 1))
rownames(prs_matrix) <- test_ids

# Process each LD threshold
for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  cat("\nProcessing LD threshold r² =", r2_str, "\n")
  
  # Process each chromosome
  for (chr in 1:22) {
    cat("  Chromosome", chr, "... ")
    
    # Construct file paths
    clump_file <- paste0(clump_pattern, chr, "_ld", r2_str, ".snplist")
    gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
    geno_prefix <- paste0(geno_pattern, chr, "_filtered")
    temp_prefix <- file.path(temp_dir, sprintf("temp_rep%s_chr%d_r2%s", rep, chr, r2_str))
    
    # Get clumped SNPs
    clumped_snps <- fread(clump_file, header = FALSE)[[1]]
    
    # Get GWAS results
    gwas <- fread(gwas_file)
        
    # Extract required samples/snps using PLINK2
    plink_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s",
      plink2_path, geno_prefix, test_ids_file, clump_file, temp_prefix
    )
    system(plink_cmd)
    
    # Read plink data
    plink_data <- read_plink(temp_prefix)
    
    # Get genotype matrix
    X <- plink_data$X # SNPs as row and individuals as column - nrow = n_snp ncol = n_sample
    snp_ids <- plink_data$bim$id
    
    # Match weights to SNPs
    weights <- gwas[ID %in% clumped_snps, .(ID, BETA)]
    matched_idx <- match(snp_ids, weights$ID)
    beta_vec <- weights$BETA[matched_idx] # n_snp * 1
    # Calculate PRS contributions
    chr_prs <-  as.vector(crossprod(X, beta_vec))   
    # Accumulate results
    sample_ids <- plink_data$fam$id
    idx <- match(sample_ids, test_ids)
    prs_matrix[idx, r2_idx] <- prs_matrix[idx, r2_idx] + chr_prs
    
    # Clean up temporary files
    unlink(c(
      paste0(temp_prefix, ".bed"),
      paste0(temp_prefix, ".bim"),
      paste0(temp_prefix, ".fam"),
      paste0(temp_prefix, ".log")
    ))
    
    cat("Done\n")
  }
}

# Save PRS results
prs_df <- data.table(IID = test_ids, prs_matrix)
output_file <- paste0(output_prefix, ".txt")
fwrite(prs_df, output_file, sep = "\t")
cat("\nPRS results saved to:", output_file, "\n")

# Clean up temporary directory
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  cat("Removed temporary directory:", temp_dir, "\n")
}