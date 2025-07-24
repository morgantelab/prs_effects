#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
bfile <- args[1]
pheno_file <- args[2]
train_ids <- args[3]
out_prefix <- args[4]
plink2_path <- args[5]

# Load data
library(data.table)
pheno <- fread(pheno_file)
train_ids_vec <- fread(train_ids, header = FALSE)[[1]]

# Prepare training phenotype
train_pheno <- pheno[IID %in% train_ids_vec, .(FID, IID, pheno)]
train_pheno_file <- tempfile(fileext = ".txt")
fwrite(train_pheno, train_pheno_file, sep = "\t", quote = FALSE)

# Run PLINK2
cmd <- sprintf(
  "%s --bfile %s --pheno %s --keep %s --glm 'no-x-sex' 'allow-no-covars' --no-input-missing-phenotype --pheno-name pheno --out %s --threads 1",
  plink2_path, bfile, train_pheno_file, train_ids, out_prefix
)
system(cmd)

# Clean up
file.remove(train_pheno_file)