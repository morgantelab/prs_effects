#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
method1 <- args[1]
method2 <- args[2]
pheno_file <- args[3]
rep <- args[4]
output1 <- args[5]
output2 <- args[6]
geno_dir <- args[7]

library(data.table)
library(dplyr)
library(BEDMatrix)

###----------------------- PRS Calculation Functions ---------------------------
calculate_prs_effect <- function(genotype, additive_effects, dominance_effects) {
  total_effect <- 0
  for (locus in names(genotype)) {
    a <- additive_effects[locus]
    d <- dominance_effects[locus]
    
    if (genotype[[locus]] == 2) {  # Homozygous alternate (e.g., AA)
      total_effect <- total_effect + a
    } else if (genotype[[locus]] == 1) {  # Heterozygous (e.g., Aa)
      total_effect <- total_effect + d
    } else if (genotype[[locus]] == 0) {  # Homozygous reference (e.g., aa)
      total_effect <- total_effect - a
    }
  }
  return(total_effect)
}

test_id <- fread(paste0("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/testing_ID_rep"
                        , rep, ".txt"))
pheno <- fread(pheno_file)
test_pheno <- pheno %>% filter(IID %in% test_id$IID)
#(1) Method 1 PRS prdiction ----------------------------------------------------

method1 <- fread("method1.txt") 

###---------------------------- PRS Calculation --------------------------------
cat("PRS Prediction using Method 1", "\n")

for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  plink_file <- paste0(geno_dir, "/ukb_chr", chr, "_filtered_200k")
  cat("PLINK file:", plink_file, "\n")
  
  # Read BIM and FAM files
  bim <- fread(paste0(plink_file, ".bim"), 
               col.names = c("chr", "SNP", "pos_cM", "pos", "alt", "ref"))
  bim <- bim[bim$SNP %in% method1$SNP[method1$chromosome == chr], ]
  
  fam <- fread(paste0(plink_file, ".fam"), 
               col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
  fam <- fam[fam$IID %in% test_id$IID, ]
  
  cat("Successfully read in bim and fam data", "\n")
  
  # Read BED file
  bed_matrix <- BEDMatrix(paste0(plink_file, ".bed"), 
                          subset = list(id = fam$IID, snp = bim$SNP))
  
  # Extract relevant genotypes
  geno <- as.matrix(bed_matrix)
  
  chr_snps <- method1 %>% filter(chromosome == chr, SNP %in% colnames(geno))
  
  # Prepare effect sizes
  additive_effects <- setNames(chr_snps$a_hat, chr_snps$SNP)
  dominance_effects <- setNames(chr_snps$d_hat, chr_snps$SNP)
  
  # Calculate PRS for each individual
  chr_prs <- apply(geno, 1, calculate_prs_effect, 
                   additive_effects = additive_effects, 
                   dominance_effects = dominance_effects)
  
  # Add to cumulative PRS
  all_prs$PRS <- all_prs$PRS + chr_prs
  
  cat("Finished processing chromosome", chr, "\n")
}

# Final PRS results
final_prs_method1 <- all_prs %>%
  filter(IID %in% test_id$IID)

# Write results to file
fwrite(final_prs_method1, output1, sep="\t")

cat("PRS calculation complete. Results written to", output1, "\n")




#(2) Method 2 PRS prdiction ----------------------------------------------------
method2 <- fread("method2.txt") 

cat("PRS Prediction using Method 2", "\n")

for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  plink_file <- paste0(geno_dir, "/ukb_chr", chr, "_filtered_200k")
  cat("PLINK file:", plink_file, "\n")
  
  # Read BIM and FAM files
  bim <- fread(paste0(plink_file, ".bim"), 
               col.names = c("chr", "SNP", "pos_cM", "pos", "alt", "ref"))
  bim <- bim[bim$SNP %in% method2$SNP[method2$chromosome == chr], ]
  
  fam <- fread(paste0(plink_file, ".fam"), 
               col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))
  fam <- fam[fam$IID %in% test_id$IID, ]
  
  cat("Successfully read in bim and fam data", "\n")
  
  # Read BED file
  bed_matrix <- BEDMatrix(paste0(plink_file, ".bed"), 
                          subset = list(id = fam$IID, snp = bim$SNP))
  
  # Extract relevant genotypes
  geno <- as.matrix(bed_matrix)
  
  chr_snps <- method2 %>% filter(chromosome == chr, SNP %in% colnames(geno))
  
  # Prepare effect sizes
  additive_effects <- setNames(chr_snps$a_hat, chr_snps$SNP)
  dominance_effects <- setNames(chr_snps$d_hat, chr_snps$SNP)
  
  # Calculate PRS for each individual
  chr_prs <- apply(geno, 1, calculate_prs_effect, 
                   additive_effects = additive_effects, 
                   dominance_effects = dominance_effects)
  
  # Add to cumulative PRS
  all_prs$PRS <- all_prs$PRS + chr_prs
  
  cat("Finished processing chromosome", chr, "\n")
}

# Final PRS results
final_prs_method2 <- all_prs %>%
  filter(IID %in% test_id$IID)

# Write results to file
fwrite(final_prs_method2, output2, sep="\t")

cat("PRS calculation complete. Results written to", output2, "\n")


