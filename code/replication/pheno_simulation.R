#!/usr/bin/env Rscript
# Get command-line arguments
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
rep <- as.integer(args[1])
seed_val <- as.integer(args[2])
geno_dir <- args[3]
results_dir <- args[4]

# Set seed for reproducibility
set.seed(seed_val)

library(data.table)
library(dplyr)
library(BEDMatrix)

h2 <- 0.5
p <- 1000
n_snps <- 1000

plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

# (1) Select random causal SNPs genome-wide --------------------------------
# Combine frequencies across chromosomes
all_freq <- rbindlist(lapply(1:22, function(chr) {
  fread(paste0(geno_dir, "/ukb_chr", chr, "_filtered_200k.afreq"),
        select = c(1, 2, 3, 4, 6),  # CHR, SNP ID, REF, ALT, ALT_FREQ
        col.names = c("chr", "snp", "ref", "alt", "alt_freq"))
}))
# Randomly select SNPs
set.seed(seed_val)
causal_snps <- all_freq[sample(.N, n_snps)][, index := .I]
causal_snps_list <- causal_snps %>% select(snp)
# Write causal SNPs to extract
write.table(causal_snps_list, 
            file = paste0(results_dir, "/causal_snps", rep, ".txt"),
            col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
# (2) Simulate effects -----------------------------------------------------
set.seed(seed_val)
a <- rnorm(p, mean = 0, sd = sqrt(2))
set.seed(seed_val)
dom_type <- rep(c("complete", "incomplete", "additive"), times = c(500, 300, 200))

d <- numeric(p)

for (i in 1:p) {
  if (dom_type[i] == "complete") {
    d[i] <- sample(c(1, -1), 1) * a[i]
  } else if (dom_type[i] == "incomplete") {
    d[i] <- runif(1, 0.2, 0.8) * sample(c(1, -1), 1) * a[i]
  } else if (dom_type[i] == "additive"){
    d[i] <- 0
  }
}

# (3) Process chromosomes ---------------------------------------------------
total_geno <- NULL
chunk_size <- 1000  # Number of samples to process at once

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr", CHR, "_filtered_200k"),
    "--extract", paste0(results_dir, "/causal_snps_rep", rep, ".txt"),
    "--make-bed",
    "--out", paste0(base_dir, "/temp_chr", CHR),
    "--memory 8000",
    "--threads 8"
  ))
  
  plink_file <- paste0(base_dir, "/temp_chr", CHR)
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
  
  # Load genotype data with BEDMatrix
  bed_path <- paste0(plink_file, ".bed")
  geno <- BEDMatrix(bed_path, simple_names = TRUE)
  
  bim_file <- paste0(plink_file, ".bim")
  bim <- fread(bim_file, col.names = c("chr", "snp", "pos_cM", "pos", "alt", "ref"))
  
  geno_matrix <- t(geno)
  colnames(geno_matrix) <- bim$snp
  
  chr_snps <- causal_snps[chr == CHR, snp]
  geno_matrix <- geno_matrix[, chr_snps, drop = FALSE]
  
  a_vals <- a[causal_snps[chr == CHR, index]]
  d_vals <- d[causal_snps[chr == CHR, index]]
  
  for (j in 1:ncol(geno_matrix)) {
    for (i in 1:nrow(geno_matrix)) {
      genotype <- geno_matrix[i, j]
      if (is.na(genotype)) next
      
      if (geno_matrix[i, j] == 0) {
        total_geno$score[i] <- total_geno$score[i] - a_vals[j]
      } else if (geno_matrix[i, j] == 1) {
        total_geno$score[i] <- total_geno$score[i] + d_vals[j]
      } else if (geno_matrix[i, j] == 2) {
        total_geno$score[i] <- total_geno$score[i] + a_vals[j]
      }
    }
  }
  
  unlink(paste0(base_dir, "/temp_chr", CHR, c(".bed", ".bim", ".fam", ".log")))
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
            paste0(results_dir, "/simulated_pheno_rep", rep, ".txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(total_geno[, c("FID", "IID", "score")],
            paste0(results_dir, "/simulated_geno_rep", rep, ".txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nSuccess! Simulated pheno and geno saved for replicate", rep, "in:", results_dir, "\n")