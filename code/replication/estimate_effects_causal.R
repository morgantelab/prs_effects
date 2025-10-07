#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
geno_pattern <- args[2] 
chunk_size <- as.integer(args[3])  # (not used now, but kept for compatibility)
output_effects <- args[4]
output_train <- args[5]
output_test <- args[6]
pheno_file <- args[7]
causal_snps_file <- args[8]        # single column: SNP IDs
plink2_path <- args[9]             # plink2 executable path
seed_val <- as.integer(args[10])
num_train <- as.integer(args[11])
num_sample <- as.integer(args[12])

cat("\nParsed arguments:\n")
cat("Replicate:", rep, "\n")
cat("Geno pattern:", geno_pattern, "\n")
cat("Chunk size:", chunk_size, "\n")
cat("Output effects:", output_effects, "\n")
cat("Output train:", output_train, "\n")
cat("Output test:", output_test, "\n")
cat("Phenotype file:", pheno_file, "\n")
cat("Causal SNPs file:", causal_snps_file, "\n")
cat("Plink2 path:", plink2_path, "\n")
cat("Seed value:", seed_val, "\n")
cat("Training size:", num_train, "\n")
cat("Sample size:", num_sample, "\n")

# Load required packages
library(data.table)
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
    
    # Additive effect
    y_AA <- pheno[AA]
    y_aa <- pheno[aa]
    mean_AA <- mean(y_AA, na.rm = TRUE)
    mean_aa <- mean(y_aa, na.rm = TRUE)
    var_AA <- var(y_AA, na.rm = TRUE)
    var_aa <- var(y_aa, na.rm = TRUE)
    
    a_hat <- (mean_AA - mean_aa) / 2
    var_a <- 0.25 * (var_AA / n_AA + var_aa / n_aa)
    df_a <- (var_a * 4)^2 / ((var_AA / n_AA)^2 / (n_AA - 1) + (var_aa / n_aa)^2 / (n_aa - 1))
    t_a <- a_hat / sqrt(var_a)
    p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)
    
    # Dominance effect
    y_Aa <- pheno[Aa]
    mean_Aa <- mean(y_Aa, na.rm = TRUE)
    var_Aa <- var(y_Aa, na.rm = TRUE)
    
    d_hat <- mean_Aa - (mean_AA + mean_aa) / 2
    var_d <- var_Aa / n_Aa + 0.25 * (var_AA / n_AA + var_aa / n_aa)
    df_d <- var_d^2 / (
      (var_Aa / n_Aa)^2 / (n_Aa - 1) +
        (var_AA / (4 * n_AA))^2 / (n_AA - 1) +
        (var_aa / (4 * n_aa))^2 / (n_aa - 1)
    )
    t_d <- d_hat / sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j, ] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}

# Read phenotype and causal SNPs
pheno <- fread(pheno_file, colClasses = c("FID" = "character", "IID" = "character"))
causal_snps <- fread(causal_snps_file, col.names = c("SNP"), header = FALSE)$SNP  # Read as single column

cat("Successfully read phenotype and", length(causal_snps), "causal SNPs for replicate", rep, "\n")

###----------------------- Split train and test set --------------------------
set.seed(seed_val)
train_ids <- sample(pheno$IID, num_train)
test_ids <- setdiff(pheno$IID, train_ids)

fwrite(pheno[IID %in% train_ids, .(FID, IID)], output_train, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(pheno[IID %in% test_ids, .(FID, IID)], output_test, sep = "\t", quote = FALSE, col.names = FALSE)

cat("Training/test IDs saved in", "\n",
    output_train, "\n", output_test, "\n")


###----------------------- Process by Chromosome --------------------------
all_results <- list()  # Collect all results

for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  # Temporary files
  temp_prefix <- paste0(tempdir(), "/chr", chr, "_causal")
  temp_bed <- paste0(temp_prefix, ".bed")
  temp_bim <- paste0(temp_prefix, ".bim")
  temp_fam <- paste0(temp_prefix, ".fam")
  snp_list <- paste0(tempdir(), "/causal_snps_chr", chr, ".txt")
  
  # Full plink prefix
  full_plink <- paste0(geno_pattern, chr, "_filtered")
  
  # Check if .bim exists
  if (!file.exists(paste0(full_plink, ".bim"))) {
    cat("Warning: No .bim file for chromosome", chr, "- skipping\n")
    next
  }
  
  # Read .bim to get SNPs
  bim <- fread(paste0(full_plink, ".bim"), 
               col.names = c("chr", "SNP", "pos_cM", "pos", "alt", "ref"))
  
  # Find causal SNPs on this chromosome
  shared_snps <- intersect(causal_snps, bim$SNP)
  if (length(shared_snps) == 0) {
    cat("No causal SNPs found on chromosome", chr, "- skipping\n")
    next
  }
  
  # Write SNP list for Plink
  fwrite(data.table(SNP = shared_snps), snp_list, col.names = FALSE, quote = FALSE)
  
  # Run Plink to extract
  plink_cmd <- sprintf(
    "%s --bfile %s --extract %s --make-bed --out %s --threads %d",
    plink2_path,
    full_plink,
    snp_list,
    temp_prefix,
    8
  )
  cat("Running:", plink_cmd, "\n")
  system(plink_cmd, intern = FALSE)
  
  if (!file.exists(temp_bed)) {
    cat("Warning: Plink2 failed for chromosome", chr, "- skipping\n")
    next
  }
  
  # Load FAM and BED
  fam <- fread(temp_fam, col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
  train_pheno <- pheno[IID %in% train_ids, ]
  pheno_vector <- train_pheno$pheno[match(fam$IID, train_pheno$IID)]
  
  bed <- BEDMatrix(temp_bed, simple_names = TRUE)
  snp_data <- bed[train_ids, , drop = FALSE]
  
  if (ncol(snp_data) == 0) {
    cat("No genotype data loaded for chr", chr, "\n")
    next
  }
  
  # Create data frame
  df <- data.frame(
    pheno = pheno_vector[match(rownames(snp_data), fam$IID)],
    as.data.frame(snp_data)
  )
  
  # Calculate effects
  results <- calculate_welch(df, colnames(snp_data))
  results_dt <- data.table(results)
  
  # Merge with bim for SNP info
  bim_sub <- bim[SNP %in% colnames(snp_data)]
  results_final <- merge(results_dt, bim_sub, by = "SNP")
  results_final <- results_final[, .(
    chromosome = chr,
    position = pos,
    SNP = SNP,
    ref_allele = ref,
    alt_allele = alt,
    a_hat, d_hat, p_a, p_d
  )]
  
  # Store in list
  all_results[[length(all_results) + 1]] <- results_final
  
  # Clean up
  unlink(paste0(temp_prefix, ".*"))
  unlink(snp_list)
  
  cat("Finished chromosome", chr, "with", nrow(results_final), "SNPs.\n")
}

###----------------------- Write All Results at Once --------------------------

if (length(all_results) > 0) {
  final_dt <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  fwrite(final_dt, output_effects, sep = "\t")
  cat("✅ Final effects written to:", output_effects, "\n")
} else {
  cat("⚠️ No results to write for replicate", rep, "\n")
}