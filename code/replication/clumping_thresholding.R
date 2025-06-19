#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
effects_file <- args[1]
output1 <- args[2]
output2 <- args[3]
rep <- args[3]
geno_dir <- args[4]

library(data.table)
library(dplyr)

plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

effects <- fread(effects_file)

training_id <- paste0("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/training_ID_rep"
                      , rep, ".txt")

# Prepare additive p-value file
add_pfile <- effects[, .(SNP, P = p_a)]
fwrite(add_pfile, "add_pfile.txt", sep = "\t", quote = FALSE)

# Prepare dominance p-value file
dom_pfile <- effects[, .(SNP, P = p_d)]
fwrite(dom_pfile, "dom_pfile.txt", sep = "\t", quote = FALSE)

# (1) Clumping and thresholding for both additive and dominance ----------------
cat("Method 1: Clumping and thresholding for both additive and dominance", "\n")

# Run clumping for additive effects
cat("Clumping and thresholding for additive effects", "\n")

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr", CHR, "_filtered_200k"),
    "--keep", paste0(training_id),
    "--clump add_pfile.txt",
    "--clump-p1", 10^(-5),
    "--clump-r2", 0.1,
    "--clump-kb", 1000,
    "--out add_clump",
    "--threads 8"
  ))
} 


# Run clumping for dominance effects
cat("Clumping and thresholding for dominance effects", "\n")

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr", CHR, "_filtered_200k"),
    "--keep", paste0(training_id),
    "--clump dom_pfile.txt",
    "--clump-p1", 10^(-5),
    "--clump-r2", 0.1,
    "--clump-kb", 1000,
    "--out add_clump",
    "--threads 8"
  ))
} 

# Read clump results
add_clump <- fread("add_clump.clumped", select = c("SNP", "CHR", "BP", "P"))
dom_clump <- fread("dom_clump.clumped", select = c("SNP", "CHR", "BP", "P"))

cat("Filtering SNPs", "\n")

all_clumped <- rbind(
  add_clump %>% mutate(model = "additive"),
  dom_clump %>% mutate(model = "dominance")
)

# Filter and combine results
method1_snps <- effects %>%
  inner_join(all_clumped, by = c("SNP", "chromosome" = "CHR", "position" = "BP")) %>%
  mutate(
    is_significant = P < p_threshold,
    keep = is_significant
  ) %>%
  filter(keep) %>%
  select(SNP, chromosome, position, a_hat, d_hat, p_a, p_d, model) %>%
  mutate(method = "method1")
fwrite(method1_snps, output1, sep = "\t")
cat("Saved", nrow(method1_snps), "SNPs for Method 1 to", output1, "\n")

# Remove duplicates (in case a SNP is significant in both models)
method1_snps <- method1_snps %>% distinct(SNP, .keep_all = TRUE)

# (2) Clumping and thresholding for dominance only -----------------------------
# Run clumping on dominance effects
cat("Method 2: Clumping and thresholding for dominance only", "\n")
dom_clump <- fread("dom_clump.clumped", select = c("SNP", "CHR", "BP", "P"))
add_sig <- effects[
  !SNP %in% dom_clump$SNP & p_a < 1e-5,
  .(SNP, chromosome, position, a_hat, d_hat, p_a, p_d)
]

method2_snps <- rbind(
  effects[SNP %in% dom_clump$SNP, 
          .(SNP, chromosome, position, a_hat, d_hat, p_a, p_d)],
  add_sig
) %>% 
  mutate(method = "method2")

fwrite(method2_snps, output2, sep = "\t")
cat("Saved", nrow(method2_snps), "SNPs for Method 2 to", output2, "\n")

