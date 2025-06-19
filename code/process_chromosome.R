#!/usr/bin/env Rscript

# Get Snakemake input and output paths
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
pheno_file <- args[1]
output_file <- args[2]
chr <- args[4]  # Chromosome number from wildcard

# Set fixed paths
base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects"
results_dir <- file.path(base_dir, "results")
# Load required packages

library(data.table)
library(dplyr)
library(genio)


# Read phenotype data
cat("Read in Phenotype files", "\n")

pheno <- fread(pheno_file, colClasses = c("character", "character", "numeric"))

# Welch's t-test function
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
        (var_aa/(4*n_aa))^2/(n_aa-1)
    )
    t_d <- d_hat/sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}
cat("Processing chromosome:", chr, "\n")

# Process chromosome data
plink_file <- paste0(base_dir, "/ukb_chr", chr, "_filtered")
cat("PLINK file:", plink_file, "\n")

plink_data <- read_plink(plink_file)
cat("Successfully read in plink data", "\n")
cat("Loaded", nrow(plink_data$fam), "samples and", nrow(plink_data$bim), "SNPs \n")

geno_df <- as.data.frame(t(plink_data$X))
cat("Successfully transpose Geno", "\n")
colnames(geno_df) <- plink_data$bim$id
geno_df <- cbind(plink_data$fam[c("fam", "id")], geno_df)
names(geno_df)[1:2] <- c("FID", "IID")
snp_info <- plink_data$bim %>%
  dplyr::select(id, chr, pos, ref, alt) %>%
  dplyr::rename(SNP = id, 
                chromosome = chr,
                position = pos,
                ref_allele = ref,
                alt_allele = alt)
cat("Successfully Extract Geno Information", "\n")

# Merge phenotype and genotype data
full_data <- merge(pheno, geno_df, by = c("FID", "IID"))
cat("Pheno and geno merger", "\n")
snp_cols <- setdiff(colnames(full_data), c("FID", "IID", "pheno"))

# Calculate effects
chr_results <- calculate_welch(full_data, snp_cols)
cat("Successfully calculated effects", "\n")

chr_results$chr <- chr

# Add SNP information
chr_results <- chr_results %>%
  left_join(snp_info, by = "SNP") %>%
  dplyr::select(chromosome, position, SNP, ref_allele, alt_allele, everything())

# Write output
cat("Writing results to:", output_file, "\n")

write.table(chr_results, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Completed processing chromosome", chr, "\n")