#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
geno_pattern <- args[2] 
chunk_size <- as.integer(args[3])
output_effects <- args[4]
output_train <- args[5]
output_test <- args[6]
pheno_file <- args[7]
seed_val <- as.integer(args[8])
num_train <- as.integer(args[9])
num_sample <- as.integer(args[10])

cat("\nConverted values:\n")
cat("Geno file pattern:", geno_pattern, "\n")
cat("rep:", rep, "\n")
cat("seed_val:", seed_val, "\n")
cat("chunk_size:", chunk_size, "\n")
cat("output_effects:", output_effects, "\n")
cat("output_train:", output_train, "\n")
cat("output_test:", output_test, "\n")
cat("pheno_file:", pheno_file, "\n")
cat("Number of Training set:", num_train, "\n")
cat("Number of Total Sample size:", num_sample, "\n")


# Load required packages
library(data.table)
library(dplyr)
library(genio)
library(BEDMatrix)

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

pheno <- fread(pheno_file, 
              colClasses = c("FID" = "character", "IID" = "character", "pheno" = "numeric"))

cat("Successfully read phenotype files for replicate", rep, "\n")

###----------------------- Split train and test set --------------------------
set.seed(seed_val)
sample_ids <- pheno$IID
train_ids <- sample(sample_ids, num_train)
test_ids <- setdiff(sample_ids, train_ids)

# Create training dataset
train_pheno <- pheno[IID %in% train_ids]
test_pheno <- pheno[IID %in% test_ids]

# Write ID files
fwrite(train_pheno[, .(FID, IID)], output_train, 
       col.names = FALSE, sep = "\t", quote = FALSE)
fwrite(test_pheno[, .(FID, IID)], output_test, 
       col.names = FALSE, sep = "\t", quote = FALSE)

cat("Created training/test ID files:\n",
    "- Training:", output_train, "\n",
    "- Test:", output_test, "\n")

###--------------- Process chromosomes using training set ----------------
# Initialize output file with header
header <- data.table(
  chromosome = integer(),
  position = integer(),
  SNP = character(),
  ref_allele = character(),
  alt_allele = character(),
  a_hat = numeric(),
  d_hat = numeric(),
  p_a = numeric(),
  p_d = numeric()
)
fwrite(header, output_effects, sep = "\t")

for (chr in 1:22) {
  cat("\n========================================")
  cat("\nProcessing chromosome", chr)
  cat("\n========================================\n")
  
  plink_file <- paste0(geno_pattern, chr, "_filtered" )
  bim_file <- paste0(plink_file, ".bim")
  fam_file <- paste0(plink_file, ".fam")
  
  bim <- fread(bim_file, col.names = c("chr", "SNP", "pos_cM", "pos", "alt", "ref"))
  fam <- fread(fam_file, col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
  
  # Match phenotype data with genotype data
  train_indices <- which(fam$IID %in% train_ids)
  pheno_vector <- train_pheno$pheno[match(fam$IID[train_indices], train_pheno$IID)]
  
  # Initialize BEDMatrix
  bed_path <- paste0(plink_file, ".bed")
  bed_matrix <- BEDMatrix(bed_path, simple_names = TRUE)
  
  # Process in chunks
  n_snps <- nrow(bim)
  n_chunks <- ceiling(n_snps / chunk_size)
  
  for (i in 1:n_chunks) {
    start_idx <- (i-1)*chunk_size + 1
    end_idx <- min(i*chunk_size, n_snps)
    snp_names <- bim$SNP[start_idx:end_idx]
    
    cat("Processing SNPs", start_idx, "to", end_idx, "\n")
    
    # Extract genotype chunk
    geno_chunk <- bed_matrix[train_indices, start_idx:end_idx, drop = FALSE]
    colnames(geno_chunk) <- snp_names
    
    # Create data chunk
    temp_df <- data.frame(pheno = pheno_vector, as.data.frame(geno_chunk))
    
    # Calculate effects
    chunk_results <- calculate_welch(temp_df, snp_names)
    
    # Prepare output
    output_chunk <- data.table(
      chromosome = bim$chr[start_idx:end_idx],
      position = bim$pos[start_idx:end_idx],
      SNP = bim$SNP[start_idx:end_idx],
      ref_allele = bim$ref[start_idx:end_idx],
      alt_allele = bim$alt[start_idx:end_idx],
      a_hat = chunk_results$a_hat,
      d_hat = chunk_results$d_hat,
      p_a = chunk_results$p_a,
      p_d = chunk_results$p_d
    )
    
    # Append to output file
    fwrite(output_chunk, output_effects, sep = "\t", append = TRUE)
  }
  cat("Finished processing chromosome", chr, "\n")
}

cat("Output file:", output_effects, "\n")