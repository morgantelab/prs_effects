#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
out_dir <- args[2]  # p-threshold specific directory
ld_thresholds_str <- args[3]
test_ids_file <- args[4]
effects_file <- args[5]
geno_pattern <- args[6]
plink_path <- args[7]
p_threshold <- args[8]  # p-threshold wildcard

cat("\nStarting PRS calculation for replicate", rep, "\n")
cat("Output directory:", out_dir, "\n")
cat("LD thresholds:", ld_thresholds_str, "\n")
cat("Test ID file:", test_ids_file, "\n")
cat("Effect file:", effects_file, "\n")
cat("Geno pattern:", geno_pattern, "\n")
cat("PLINK path:", plink_path, "\n")
cat("P-threshold:", p_threshold, "\n")

library(data.table)
library(genio)

# Create temp directory
temp_dir <- file.path(out_dir, sprintf("temp_rep%d", rep))
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created temp directory:", temp_dir, "\n")
}

# Convert LD thresholds to numeric vector
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))
cat("LD thresholds:", ld_thresholds, "\n")

# Load test set IDs
test_ids <- fread(test_ids_file, header = FALSE)[[1]]

# Load effect estimates
eff <- fread(effects_file)

# Prepare PRS storage
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = length(ld_thresholds))
colnames(prs_matrix) <- sprintf("r2_%s", format(ld_thresholds, nsmall = 1))

# Process each LD threshold
for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  cat("  LD threshold r² =", r2_str, "\n")
  
  # Form clump file path
  clump_file <- file.path(out_dir, sprintf("clumped_snps_rep%d_r2_%.1f.txt", rep, r2))
  cat("  Clump file:", clump_file, "\n")
  
  # Load clumped SNPs
  if (!file.exists(clump_file)) {
    cat("  Clump file not found:", clump_file, "\n")
    next
  }
  
  clumped <- fread(clump_file)
  
  # Handle column naming differences
  if ("#CHROM" %in% names(clumped)) setnames(clumped, "#CHROM", "CHR")
  if ("ID" %in% names(clumped)) setnames(clumped, "ID", "SNP")
  if("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")
  
  # Get effect sizes for clumped SNPs
  clumped_effects <- merge(clumped[, .(SNP, CHR)], 
                           eff, 
                           by = c("SNP", "CHR"))
  
  # Process each chromosome
  for (chr in 1:22) {
    cat("    Processing chromosome", chr, "\n")
    
    # Get SNPs for this chromosome
    chr_snps <- clumped_effects[CHR == chr, 
                                .(SNP, ref_allele, alt_allele, a_hat, d_hat)]
    
    if (nrow(chr_snps) == 0) {
      cat("      No SNPs for chromosome", chr, "\n")
      next
    }
    
    # Write SNP list
    snp_list_file <- file.path(temp_dir, sprintf("chr%d_snps.txt", chr))
    fwrite(chr_snps[, .(SNP)], snp_list_file, col.names = FALSE)  
    
    # Prepare genotype files
    plink_file <- paste0(geno_pattern, chr, "_filtered")
    temp_prefix <- file.path(temp_dir, sprintf("chr%d_r2_%s", chr, r2_str))
    
    plink_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s --memory 2000 --threads 40",
      plink_path,
      plink_file,
      test_ids_file,
      snp_list_file,
      temp_prefix
    )
    system(plink_cmd)
    
    # Read genotype data
    plink_data <- genio::read_plink(temp_prefix, verbose = FALSE)
    fam_ids <- plink_data$fam$id
    snp_ids <- plink_data$bim$id
    X <- t(plink_data$X)  
    
    # Find test individuals and SNPs
    test_logical <- fam_ids %in% test_ids
    snp_logical <- snp_ids %in% chr_snps$SNP
    
    # Subset genotype matrix
    X_sub <- X[test_logical, snp_logical, drop = FALSE]
    matched_iids <- fam_ids[test_logical]
    matched_snp_ids <- snp_ids[snp_logical]
    
    # Get SNP effects for the subset
    snp_effects <- chr_snps[match(matched_snp_ids, chr_snps$SNP), ]
    
    # Add PRS contributions
    for (j in 1:ncol(X_sub)) {
      g <- X_sub[, j]
      a <- snp_effects$a_hat[j]
      d <- snp_effects$d_hat[j]
      
      # Calculate contributions
      contrib <- rep(0, length(g))
      contrib[g == 0] <- a  # Homozygous reference
      contrib[g == 1] <- d   # Heterozygous
      contrib[g == 2] <- -a   # Homozygous alternate
      
      # Find positions in overall matrix
      row_indices <- match(matched_iids, test_ids)
      
      # Accumulate to PRS matrix
      prs_matrix[row_indices, r2_idx] <- prs_matrix[row_indices, r2_idx] + contrib
    }
  }
}

# Cleanup temporary files
unlink(temp_dir, recursive = TRUE)

# Save PRS results
prs_df <- data.table(iid = test_ids, prs_matrix)
output_file <- file.path(out_dir, sprintf("test_prs_scan_rep%d.txt", rep))
fwrite(prs_df, output_file, sep = "\t")
cat("Saved PRS results to:", output_file, "\n")