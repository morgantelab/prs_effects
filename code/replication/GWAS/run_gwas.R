#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
rep <- args[1]
results_dir <- args[2]
ld_thresholds <- args[3]
train_ids <- args[4]
pheno_file <- args[5]
out_prefix <- args[6]
plink2_path <- args[7]
geno_pattern <- args[8]

cat("\nInput Parametes:\n")
cat("Geno file pattern:", geno_pattern, "\n")
cat("rep:", rep, "\n")
cat("ld_thresholds:", ld_thresholds, "\n")
cat("Training Files:", train_ids, "\n")
cat("Pheno Files:", pheno_file, "\n")
cat("Output Prefix:", out_prefix, "\n")
cat("Results directory:", results_dir, "\n")

# Load data
library(data.table)
pheno <- fread(pheno_file)
train_ids_vec <- fread(train_ids, header = FALSE)[[1]]

# Prepare training phenotype
train_pheno <- pheno[IID %in% train_ids_vec, .(FID, IID, pheno)]
train_pheno_file <- paste0(results_dir, "train_pheno_rep", rep, ".txt")
fwrite(train_pheno, train_pheno_file, sep = "\t", quote = FALSE)

# Run PLINK2
chromosomes <- 1:22

for(chr in chromosomes) {
  bfile <- paste0(geno_pattern, chr, "_filtered")
  output_file <- paste0(out_prefix, chr)
  
  plink2_cmd <- sprintf(
    "%s --bfile %s --pheno %s --keep %s --glm no-x-sex allow-no-covars --no-input-missing-phenotype --pheno-name pheno --out %s --threads 1",
    plink2_path, bfile, train_pheno_file, train_ids, output_file
  )
  
  cat("\nRunning PLINK2 command for chromosome", chr, ":\n", plink2_cmd, "\n")
  system(plink2_cmd)
}

# Clean up
file.remove(train_pheno_file)