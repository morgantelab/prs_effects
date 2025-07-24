#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
results_dir <- args[1]
geno_dir <- args[2]
gwas_dir <- args[3]
temp_dir <- args[4]
plink2_path <- args[5]
ld_thresholds <- as.numeric(unlist(strsplit(args[6], ",")))
rep <- args[7]
test_ids_file <- args[8]

library(data.table)
library(genio)

# Load test IDs
test_ids <- fread(test_ids_file, header = FALSE)[[1]]

# Initialize PRS matrix
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = length(ld_thresholds))
colnames(prs_matrix) <- sprintf("r2_%s", format(ld_thresholds, nsmall = 1))
rownames(prs_matrix) <- test_ids

# Process each LD threshold
for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  
  # Process each chromosome
  for (chr in 1:22) {
    # Get clumped SNPs
    clump_file <- file.path(gwas_dir, sprintf("gwas_clump_rep%s_chr%d_ld%s.snplist", rep, chr, r2_str))
    clumped_snps <- fread(clump_file, header = FALSE)[[1]]
    
    # Get GWAS results
    gwas_file <- file.path(gwas_dir, sprintf("gwas_rep%s_chr%d.pheno.glm.linear", rep, chr))
    gwas <- fread(gwas_file)
    
    # Prepare genotype data
    geno_prefix <- file.path(geno_dir, sprintf("ukb_chr%d_filtered", chr))
    temp_prefix <- file.path(temp_dir, sprintf("temp_rep%s_chr%d_ld%s", rep, chr, r2_str))
    
    # Extract required samples/snps
    cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s",
      plink2_path, geno_prefix, test_ids_file, clump_file, temp_prefix
    )
    system(cmd)
    
    # Read plink data
    plink_data <- read_plink(temp_prefix)
    
    # Get genotype matrix
    X <- plink_data$X
    X[is.na(X)] <- 0
    snp_ids <- plink_data$bim$id
    
    # Match weights to SNPs
    weights <- gwas[ID %in% clumped_snps, .(ID, BETA)]
    matched_idx <- match(snp_ids, weights$ID)
    beta_vec <- weights$BETA[matched_idx]
    
    # Calculate PRS contributions
    chr_prs <- as.vector(X %*% beta_vec)
    
    # Accumulate results
    sample_ids <- plink_data$fam$id
    idx <- match(sample_ids, test_ids)
    prs_matrix[idx, r2_idx] <- prs_matrix[idx, r2_idx] + chr_prs
    
    # Clean up temporary files
    unlink(c(
      paste0(temp_prefix, ".bed"),
      paste0(temp_prefix, ".bim"),
      paste0(temp_prefix, ".fam")
    ))
  }
}

# Save PRS results
prs_df <- data.table(IID = test_ids, prs_matrix)
output_file <- file.path(gwas_dir, sprintf("prs_results_beta_rep%s.txt", rep))
fwrite(prs_df, output_file, sep = "\t")