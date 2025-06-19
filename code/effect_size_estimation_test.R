#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(genio)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chr <- as.integer(args[1])

# Load configuration
config <- yaml::read_yaml("config.yaml")
base_dir <- config$base_dir
results_dir <- config$results_dir
chr_dir <- file.path(results_dir, "chr_effects_test")
if (!dir.exists(chr_dir)) dir.create(chr_dir, recursive = TRUE)

###----------------------- Effect Calculation Function ------------------------
calculate_welch <- function(data, snp_cols) {
  n_snps <- length(snp_cols)
  results <- data.frame(
    SNP = snp_cols,
    a_hat = numeric(n_snps),
    d_hat = numeric(n_snps),
    p_a = numeric(n_snps),
    p_d = numeric(n_snps)
  )
  
  for(j in 1:n_snps) {
    snp <- snp_cols[j]
    X <- data[[snp]]
    pheno <- data$pheno
    
    AA <- X == 2
    Aa <- X == 1
    aa <- X == 0
    
    n_AA <- sum(AA, na.rm = TRUE)
    n_Aa <- sum(Aa, na.rm = TRUE)
    n_aa <- sum(aa, na.rm = TRUE)
    
    # Additive effect calculations
    y_AA <- pheno[AA]
    y_aa <- pheno[aa]
    
    mean_AA <- mean(y_AA, na.rm = TRUE)
    mean_aa <- mean(y_aa, na.rm = TRUE)
    var_AA <- var(y_AA, na.rm = TRUE)
    var_aa <- var(y_aa, na.rm = TRUE)
    
    a_hat <- (mean_AA - mean_aa)/2
    var_a <- 0.25*(var_AA/n_AA + var_aa/n_aa)
    df_a <- (var_AA/n_AA + var_aa/n_aa)^2 / 
      ((var_AA/n_AA)^2/(n_AA-1) + (var_aa/n_aa)^2/(n_aa-1))
    t_a <- a_hat/sqrt(var_a)
    p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)
    
    # Dominance effect calculations
    y_Aa <- pheno[Aa]
    mean_Aa <- mean(y_Aa, na.rm = TRUE)
    var_Aa <- var(y_Aa, na.rm = TRUE)
    
    d_hat <- mean_Aa - (mean_AA + mean_aa)/2
    var_d <- var_Aa/n_Aa + 0.25*(var_AA/n_AA + var_aa/n_aa)
    df_d <- (var_d)^2 / (
      (var_Aa/n_Aa)^2/(n_Aa-1) + 
        (var_AA/(4*n_AA))^2/(n_AA-1) + 
        (var_aa/(4*n_aa))^2/(n_aa-1))
    t_d <- d_hat/sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}

###----------------------- Process Single Chromosome --------------------------
cat("Processing chromosome", chr, "\n")

# Read phenotype data
pheno <- fread(file.path(results_dir, "simulated_pheno_large.txt"),
               colClasses = c("FID" = "character", "IID" = "character", "pheno" = "numeric"))

# Read PLINK files
plink_file <- file.path(base_dir, paste0("ukb_chr", chr, "_filtered"))
plink_data <- read_plink(plink_file)

# Create genotype data frame
geno_df <- as.data.frame(t(plink_data$X))
colnames(geno_df) <- plink_data$bim$id
geno_df <- cbind(plink_data$fam[, c("fam", "id")], geno_df)
names(geno_df)[1:2] <- c("FID", "IID")

# Extract SNP info
snp_info <- plink_data$bim %>%
  dplyr::select(id, chr, pos, ref, alt) %>%
  dplyr::rename(SNP = id, 
                chromosome = chr,
                position = pos,
                ref_allele = ref,
                alt_allele = alt)

# Merge with phenotype data
full_data <- merge(pheno, geno_df, by = c("FID", "IID"))
snp_cols <- setdiff(colnames(full_data), c("FID", "IID", "pheno"))

# Calculate Welch statistics
chr_results <- calculate_welch(full_data, snp_cols)
chr_results$chr <- chr

# Add SNP information
chr_results <- left_join(chr_results, snp_info, by = "SNP") %>%
  dplyr::select(chromosome, position, SNP, ref_allele, alt_allele, everything())

# Save results
output_file <- file.path(chr_dir, paste0("chr", chr, "_effects.txt"))
fwrite(chr_results, output_file, sep = "\t", quote = FALSE)
cat("Saved results for chromosome", chr, "\n")