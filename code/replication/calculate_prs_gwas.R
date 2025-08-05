#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments with new parameter order
geno_pattern <- args[1]        # Genotype file pattern
gwas_pattern <- args[2]        # GWAS results pattern
clump_pattern <- args[3]       # Clumping results pattern
test_ids_file <- args[4]       # File with test IDs
plink2_path <- args[5]         # Path to PLINK2
ld_thresholds_str <- args[6]   # Comma-separated LD thresholds
out_dir <- args[7]             # Output directory (p-threshold specific)
rep <- args[8]                 # Replicate number
p_threshold <- args[9]         # p-threshold wildcard
output_file <- args[10]        # Output file name

cat("\nStarting PRS calculation for replicate", rep, "p-threshold", p_threshold, "\n")
cat("Genotype pattern:", geno_pattern, "\n")
cat("GWAS pattern:", gwas_pattern, "\n")
cat("Clump pattern:", clump_pattern, "\n")
cat("Test IDs file:", test_ids_file, "\n")
cat("PLINK2 path:", plink2_path, "\n")
cat("LD thresholds:", ld_thresholds_str, "\n")
cat("Output directory:", out_dir, "\n")
cat("Output file:", output_file, "\n")

# Convert LD thresholds to numeric vector
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))
cat("LD thresholds:", ld_thresholds, "\n")

# Create temporary directory within p-threshold specific directory
temp_dir <- file.path(out_dir, paste0("temp_prs_rep", rep))
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
prs_matrix <- matrix(0, 
                     nrow = length(test_ids), 
                     ncol = length(ld_thresholds),
                     dimnames = list(test_ids, paste0("r2_", format(ld_thresholds, nsmall = 1))))
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
    
    # Construct file paths with new directory structure
    clump_file <- paste0(clump_pattern, chr, "_r2_", r2_str, ".snplist")
    gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
    geno_prefix <- paste0(geno_pattern, chr, "_filtered")
    temp_prefix <- file.path(temp_dir, sprintf("chr%d_r2%s", chr, r2_str))
    
    # Check if clump file exists
    if (!file.exists(clump_file)) {
      cat("  Clump file not found:", clump_file, "\n")
      next
    }
    
    # Get clumped SNPs
    clumped_snps <- tryCatch({
      fread(clump_file, header = FALSE)[[1]]
    }, error = function(e) {
      cat("  Error reading clump file:", clump_file, "\n")
      character(0)
    })
    
    if (length(clumped_snps) == 0) {
      cat("  No clumped SNPs for chromosome", chr, "\n")
      next
    }
    
    # Get GWAS results
    gwas <- tryCatch({
      fread(gwas_file)
    }, error = function(e) {
      cat("  Error reading GWAS file:", gwas_file, "\n")
      NULL
    })
    
    if (is.null(gwas)) next
    
    # Extract required samples/snps using PLINK2
    plink_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s",
      plink2_path, geno_prefix, test_ids_file, clump_file, temp_prefix
    )
    system(plink_cmd)
    
    # Check if PLINK output exists
    bed_file <- paste0(temp_prefix, ".bed")
    if (!file.exists(bed_file)) {
      cat("  PLINK output not found for chromosome", chr, "\n")
      next
    }
    
    # Read plink data
    plink_data <- tryCatch({
      read_plink(temp_prefix)
    }, error = function(e) {
      cat("  Error reading PLINK files for chromosome", chr, "\n")
      NULL
    })
    
    # Get genotype matrix
    X <- plink_data$X   # SNPs as rows, individuals as columns
    snp_ids <- plink_data$bim$id
    X[is.na(X)] <- 0
    # ⚠️ WARNING: Treating missing genotypes as 0 (homozygous reference)
    # This may bias PRS downward if missingness is non-random or MAF > 0.
    # Match weights to SNPs
    weights <- gwas[ID %in% clumped_snps, .(ID, BETA)]
    matched_idx <- match(snp_ids, weights$ID)
    beta_vec <- weights$BETA[matched_idx]
    
    # Calculate PRS contributions (A2 is effect allele)
    chr_prs <- as.vector(crossprod(X, beta_vec))
    
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
    
    cat("Done -", length(beta_vec), "SNPs\n")
  }
}

# Save PRS results
prs_df <- data.table(IID = test_ids, prs_matrix)
fwrite(prs_df, output_file, sep = "\t")
cat("\nPRS results saved to:", output_file, "\n")

# Clean up temporary directory
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  cat("Removed temporary directory:", temp_dir, "\n")
}