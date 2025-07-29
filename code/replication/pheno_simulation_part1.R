#!/usr/bin/env Rscript
# Get command-line arguments
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

cat("Received arguments:\n")
print(args)

bim_pattern <- args[1]  # Get the space-separated BIM file paths
rep <- as.integer(args[2])
seed_val <- as.integer(args[3])
geno_dir <- args[4]
results_dir <- args[5]
h2 <- as.numeric(args[6])
num_loci <- as.integer(args[7])
plink2_path <- args[8]
output_causal_snps <- args[9]
output_pheno <- args[10]
output_geno <- args[11]
num_additive <- as.integer(args[12])
num_incom_dom <- as.integer(args[13])
num_com_dom <- as.integer(args[14])

cat("\nConverted values:\n")
cat("BIM file pattern:", bim_pattern, "\n")
cat("rep:", rep, "\n")
cat("seed_val:", seed_val, "\n")
cat("h2:", h2, "\n")
cat("num_loci:", num_loci, "\n")
cat("num_additive:", num_additive, "\n")
cat("num_incom_dom:", num_incom_dom, "\n")
cat("num_com_dom:", num_com_dom, "\n")

library(data.table)
library(dplyr)
library(BEDMatrix)
library(genio)

# (1) Select random causal SNPs genome-wide --------------------------------
# Combine frequencies across chromosomes
bim_files <- paste0(bim_pattern, 1:22, "_filtered.bim")
cat("\nBIM files:\n")
print(bim_files)

all_snp <- unlist(lapply(bim_files, function(bim_file) {
  fread(bim_file, select = 2, header = FALSE)[[1]]
}))

# Randomly select SNPs
set.seed(seed_val)
selected_indices <- sample(length(all_snp), num_loci)

causal_snps <- data.frame(
  snp = all_snp[selected_indices])

# Write causal SNPs to extract
write.table(causal_snps, 
            file = output_causal_snps,
            col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)

# (2) Simulate effects -----------------------------------------------------
set.seed(seed_val)
a <- rnorm(num_loci, mean = 0, sd = sqrt(2))

set.seed(seed_val)
dom_type <- rep(c("complete", "incomplete", "additive"), 
               times = c(num_com_dom, num_incom_dom, num_additive))

d <- numeric(num_loci)
# For complete dominance: random sign * a
complete_idx <- which(dom_type == "complete")
d[complete_idx] <- a[complete_idx] * sample(c(-1, 1), length(complete_idx), replace = TRUE)

# For incomplete dominance: random sign * random factor * a
incomplete_idx <- which(dom_type == "incomplete")
d[incomplete_idx] <- a[incomplete_idx] * 
                     sample(c(-1, 1), length(incomplete_idx), replace = TRUE) *
                     runif(length(incomplete_idx), 0.2, 0.8)

additive_idx <- which(dom_type == "additive")
d[additive_idx] <- 0

# (3) Process chromosomes ---------------------------------------------------
total_geno <- NULL

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr", CHR, "_filtered"),
    "--extract", output_causal_snps,
    "--make-bed",
    "--out", paste0(geno_dir, "/temp_chr", CHR, "_", rep)
  ))
  
  plink_file <- paste0(geno_dir, "/temp_chr", CHR, "_", rep)
  bim_file <- paste0(plink_file, ".bim")
  fam_file <- paste0(plink_file, ".fam")
  
  # Read sample information
  fam <- fread(fam_file, col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
  
  if(is.null(total_geno)) {
    total_geno <- data.frame(
      FID = fam$FID,
      IID = fam$IID,
      score = 0
    )
  }
  
  # Load genotype data with BEDMatrix (no transpose)
  bed_path <- paste0(plink_file, ".bed")
  geno <- BEDMatrix(bed_path)
  cat("Loading Genotype Matrix \n")
  # Read BIM file
  bim <- fread(bim_file, col.names = c("chr", "snp", "pos_cM", "pos", "alt", "ref"))
  cat("Loading BIM file \n")
  
  causal_snps <- fread(output_causal_snps, 
                       header = FALSE, col.names = "snp")
  # Get effect sizes for each SNP in this chromosome
  match_idx <- match(bim$snp, causal_snps$snp)
  a_vals <- a[match_idx]
  d_vals <- d[match_idx]
  
  # Calculate genetic score using vectorized operations
  for (j in 1:ncol(geno)) {
    g <- geno[, j]  # Genotypes for all individuals for this SNP
    
    # Calculate contributions using vectorization
    contrib <- case_when(
      g == 0 ~ -a_vals[j],
      g == 1 ~ d_vals[j],
      g == 2 ~ a_vals[j],
      TRUE ~ 0
    )
    
    total_geno$score <- total_geno$score + contrib
  }
  
  unlink(paste0(geno_dir, "/temp_chr", CHR, "_", rep, c(".bed", ".bim", ".fam", ".log")))
}


# (4) Final phenotype calculation --------------------------------------------
var_g <- var(total_geno$score)
var_e <- (1/h2 - 1) * var_g
set.seed(seed_val)
total_geno$pheno <- total_geno$score + rnorm(nrow(total_geno), 0, sqrt(var_e))

pheno <- data.frame(
  FID = as.character(total_geno$FID),
  IID = as.character(total_geno$IID),
  pheno = total_geno$pheno
)

# Save results
write.table(pheno,
            output_pheno,
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(total_geno[, c("FID", "IID", "score")],
            output_geno,
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nSuccess! Simulated pheno and geno saved for replicate", rep, "in:", results_dir, "\n")