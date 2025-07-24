#!/usr/bin/env Rscript
# Get command-line arguments
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
rep <- as.integer(args[1])
seed_val <- as.integer(args[2])
geno_dir <- args[3]
results_dir <- args[4]


# Packages and paramteres-------------------------------------------------------
library(data.table)
library(dplyr)
library(BEDMatrix)
library(genio)

h2 <- 0.5
p <- 1000
n_snps <- 1000

plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
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

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr", CHR, "_filtered"),
    "--extract", paste0(results_dir, "/causal_snps", rep, ".txt"),
    "--make-bed",
    "--out", paste0(geno_dir, "/temp_chr", CHR, "_", rep),
    "--memory 8000",
    "--threads 24"
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
  
  # Get effect sizes for each SNP in this chromosome
  match_idx <- match(bim$snp, causal_snps$snp)
  a_vals <- a[match_idx]
  d_vals <- d[match_idx]
  
  # Calculate genetic score using vectorized operations
  for (j in 1:length(matched_snps)) {
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
            paste0(results_dir, "/simulated_pheno_rep", rep, ".txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(total_geno[, c("FID", "IID", "score")],
            paste0(results_dir, "/simulated_geno_rep", rep, ".txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nSuccess! Simulated pheno and geno saved for replicate", rep, "in:", results_dir, "\n")