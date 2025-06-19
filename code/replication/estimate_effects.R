#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(genio)
library(BEDMatrix)
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
pheno_file <- args[1]
output_file <- args[2]
chr <- as.numeric(args[3])
geno_dir <- args[4]
rep <- args[5]


print(.libPaths())

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


# Read phenotype data
pheno <- fread(pheno_file, 
               colClasses = c("character", "character", "numeric"))

cat("Successfully Read in Phenotype files", "\n")

###----------------------- Split train and test set --------------------------
set.seed(1000 + as.integer(rep))
sample_ids <- pheno_full$IID
train_ids <- sample(sample_ids, 180000)
test_ids <- setdiff(sample_ids, train_ids)

# Create training dataset
train_pheno <- pheno[IID %in% train_ids]
test_pheno <- pheno[IID %in% test_ids]
trian_pheno_ID <- train_pheno %>% select(FID, IID)
write.table(train_pheno_ID, 
            file = paste0(results_dir, "/training_ID_rep", rep, ".txt"),
            col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
test_pheno_ID <- test_pheno %>% select(FID, IID)
write.table(test_pheno_ID, 
            file = paste0(results_dir, "/test_ID_rep", rep, ".txt"),
            col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)

###--------------- Process Single Chromosome using training set ----------------
# Process chromosome data
cat("Processing chromosome", chr, "\n")
plink_file <- paste0(geno_dir, "/ukb_chr", chr, "_filtered_200k")
cat("PLINK file:", plink_file, "\n")
bim_file <- paste0(plink_file, ".bim")
fam_file <- paste0(plink_file, ".fam")

bim <- fread(bim_file, col.names = c("chr", "SNP", "pos_cM", "pos", "alt", "ref"))
cat("Successfully read in bim data", "\n")

fam <- fread(fam_file, col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
cat("Successfully read in fam data", "\n")

pheno_vector <- train_pheno$pheno[match(fam$IID, train_pheno$IID)]
train_indices <- which(fam$IID %in% train_ids)

bed_path <- paste0(plink_file, ".bed")
bed_matrix <- BEDMatrix(bed_path, simple_names = TRUE)

n_snps <- nrow(bim)
chunk_size <- 1000  # Process 1000 SNPs per chunk
n_chunks <- ceiling(n_snps / chunk_size)

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

fwrite(header, output_file , sep = "\t")

for (i in 1:n_chunks) {
  start_idx <- (i-1)*chunk_size + 1
  end_idx <- min(i*chunk_size, n_snps)
  cat("Processing SNPs", start_idx, "to", end_idx, "\n")
  
  # Extract genotype chunk (all individuals)
  geno_chunk <- bed_matrix[train_indices, start_idx:end_idx, drop = FALSE]
  snp_names <- bim$SNP[start_idx:end_idx]
  
  # Create data chunk
  temp_df <- data.frame(pheno = pheno_vector)
  temp_df <- cbind(temp_df, as.data.frame(geno_chunk))
  colnames(temp_df)[-1] <- snp_names
  
  # Calculate effects
  chunk_results <- calculate_welch(temp_df, snp_names)
  
  # Prepare output with SNP info (subset rows, not columns!)
  chunk_bim <- bim[start_idx:end_idx, .(chr, SNP, pos, ref, alt)]
  output_chunk <- data.table(
    chromosome = chunk_bim$chr,
    position = chunk_bim$pos,
    SNP = chunk_bim$SNP,
    ref_allele = chunk_bim$ref,
    alt_allele = chunk_bim$alt,
    a_hat = chunk_results$a_hat,
    d_hat = chunk_results$d_hat,
    p_a = chunk_results$p_a,
    p_d = chunk_results$p_d
  )
  
  # Append to output file
  fwrite(output_chunk, output_file, sep = "\t", append = TRUE)
  
  # Explicit memory management
  rm(geno_chunk, temp_df, chunk_results, output_chunk)
  gc()
}

cat("Finished processing chromosome", chr, "for replicate", rep, "\n")
