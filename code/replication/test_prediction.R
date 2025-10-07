#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
out_dir <- args[2]              # p-threshold specific directory
ld_thresholds_str <- args[3]
test_ids_file <- args[4]
effects_file <- args[5]
geno_pattern <- args[6]
plink_path <- args[7]
p_threshold <- as.numeric(args[8])  # p-threshold wildcard

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

# Load test set IDs (first column assumed to be IID if single-col file)
test_ids <- fread(test_ids_file, header = FALSE)[[1]]
cat("✅ Loaded", length(test_ids), "test individuals\n")

# Load effect estimates
eff <- fread(effects_file)

# Harmonize column names possibly coming from different sources
if ("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")
if ("#CHROM" %in% names(eff)) setnames(eff, "#CHROM", "CHR")
if (!("CHR" %in% names(eff))) stop("Effects table must contain a 'CHR' column.")
if (!("SNP" %in% names(eff))) stop("Effects table must contain a 'SNP' column.")

# Safe addition function: keep NA but do not affect other SNPs
safe_add <- function(x, y) {
  out <- x
  both <- !is.na(x) & !is.na(y)
  only_y <-  is.na(x) & !is.na(y)
  out[both] <- x[both] + y[both]
  out[only_y] <- y[only_y]
  out
}

# Prepare PRS storage
ncol_prs <- length(ld_thresholds)
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = ncol_prs,
                     dimnames = list(test_ids, paste0("r2_", format(ld_thresholds, nsmall = 1))))

# Process each LD threshold
for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  cat("\n=== LD threshold r² =", r2_str, "===\n")
  
  # Initialize PRS accumulator for this r2 across all chromosomes
  prs_vector <- numeric(length(test_ids))
  prs_vector[] <- 0
  
  # Form clump file path
  clump_file <- file.path(out_dir, sprintf("clumped_snps_rep%d_r2_%.1f.txt", rep, r2))
  cat("Clump file:", clump_file, "\n")
  
  # Load clumped SNPs
  if (!file.exists(clump_file)) {
    cat("  ⚠️  Clump file not found, skipping r2 =", r2_str, "\n")
    next
  }
  clumped <- fread(clump_file)
  if ("#CHROM" %in% names(clumped)) setnames(clumped, "#CHROM", "CHR")
  if ("ID" %in% names(clumped)) setnames(clumped, "ID", "SNP")
  if (!("CHR" %in% names(clumped)) || !("SNP" %in% names(clumped))) {
    stop("Clump file must contain CHR and SNP columns.")
  }
  
  # Get effect sizes for clumped SNPs (limit columns we need)
  clumped_effects <- merge(
    clumped[, .(SNP, CHR)],
    eff[, .(SNP, CHR, ref_allele, alt_allele, a_hat, d_hat, p_a, p_d)],
    by = c("SNP", "CHR")
  )
  if (nrow(clumped_effects) == 0) {
    cat("  No overlapping SNPs between clumped list and effects for r2 =", r2_str, "\n")
    prs_matrix[, r2_idx] <- prs_vector
    next
  }
  
  # Process each chromosome
  for (chr in 1:22) {
    cat("  Processing chromosome", chr, "...\n")
    
    # SNPs for this chromosome
    chr_snps <- clumped_effects[CHR == chr,
                                .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_a, p_d)]
    if (nrow(chr_snps) == 0) {
      cat("    No SNPs for chromosome", chr, "\n")
      next
    }
    
    # Write SNP list for PLINK extraction
    snp_list_file <- file.path(temp_dir, sprintf("chr%d_snps.txt", chr))
    fwrite(chr_snps[, .(SNP)], snp_list_file, col.names = FALSE)
    
    # Prepare genotype files
    plink_file <- paste0(geno_pattern, chr, "_filtered")
    temp_prefix <- file.path(temp_dir, sprintf("chr%d_r2_%s", chr, r2_str))
    
    plink_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s --memory 2000 --threads 10",
      plink_path,
      plink_file,
      test_ids_file,
      snp_list_file,
      temp_prefix
    )
    system(plink_cmd)
    
    # If PLINK made no variants for this chr/r2, skip
    bed_file <- paste0(temp_prefix, ".bed")
    bim_file <- paste0(temp_prefix, ".bim")
    fam_file <- paste0(temp_prefix, ".fam")
    if (!file.exists(bed_file) || !file.exists(bim_file) || !file.exists(fam_file)) {
      cat("    No extracted genotype data for chr", chr, "\n")
      next
    }
    
    # Read genotype data
    plink_data <- genio::read_plink(temp_prefix, verbose = FALSE)
    fam_ids <- plink_data$fam$id           # individual IDs in extracted file
    snp_ids <- plink_data$bim$id           # SNP IDs in extracted file
    X <- t(plink_data$X)                   # make Individual x SNP
    
    # Ensure the order of individuals matches test_ids (it should if --keep is ordered)
    # But reindex defensively:
    if (!all(fam_ids == test_ids)) {
      # Map fam_ids to the requested test_ids order if both are the same set
      idx_map <- match(test_ids, fam_ids)
      if (any(is.na(idx_map))) {
        stop("Mismatch between test_ids and FAM IDs after --keep.")
      }
      X <- X[idx_map, , drop = FALSE]
      fam_ids <- fam_ids[idx_map]
    }
    
    # Align effects to the order of columns in X
    snp_effects <- chr_snps[match(snp_ids, chr_snps$SNP), ]
    if (any(is.na(snp_effects$SNP))) {
      # Some SNPs in the PLINK output not found in effect table; drop those columns
      keep <- !is.na(snp_effects$SNP)
      X <- X[, keep, drop = FALSE]
      snp_ids <- snp_ids[keep]
      snp_effects <- snp_effects[keep]
    }
    if (length(snp_ids) == 0) {
      cat("    After alignment, 0 SNPs remain for chr", chr, "\n")
      next
    }
    
    # Indices we’ll use
    snp_indices <- seq_along(snp_ids)
    test_indices <- seq_len(nrow(X))
    
    # Accumulate PRS contributions for this chromosome
    prs_contrib <- numeric(nrow(X))
    prs_contrib[] <- 0
    
    for (j in snp_indices) {
      g <- X[test_indices, j]  # genotypes 0/1/2/NA
      
      a <- snp_effects$a_hat[j]
      d <- snp_effects$d_hat[j]
      p_d_val <- snp_effects$p_d[j]
      
      # Contribution per genotype
      contrib <- numeric(length(g))
      # default all zero; then fill by genotype (NAs will remain 0 unless you prefer NA)
      if (!is.na(p_d_val) && p_d_val < as.numeric(p_threshold)) {
        # include dominance
        contrib[g == 0] <- -a
        contrib[g == 1] <- d
        contrib[g == 2] <- a
      } else {
        # additive only
        contrib[g == 0] <- -a
        contrib[g == 1] <- 0
        contrib[g == 2] <- a
      }
      
      # Keep missing genotypes as 0; if you prefer to propagate NA, replace with:
      # contrib[is.na(g)] <- NA_real_
      prs_contrib <- safe_add(prs_contrib, contrib)
    }
    
    # Add this chromosome’s contribution to the running PRS for this r2
    prs_vector <- safe_add(prs_vector, prs_contrib)
  } # end chr loop
  
  # Store the PRS for this r2 into the matrix (in test_ids order)
  prs_matrix[, r2_idx] <- prs_vector
} # end r2 loop

# Cleanup temporary files
unlink(temp_dir, recursive = TRUE)

# Save PRS results
prs_df <- data.table(IID = test_ids, prs_matrix)
output_file <- file.path(out_dir, sprintf("test_prs_scan_rep%d.txt", rep))
fwrite(prs_df, output_file, sep = "\t")
cat("Saved PRS results to:", output_file, "\n")
