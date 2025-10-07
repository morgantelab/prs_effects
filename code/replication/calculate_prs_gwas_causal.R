#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

geno_pattern <- args[1]
gwas_pattern <- args[2]
test_ids_file <- args[3]
plink2_path <- args[4]
outdir <- args[5]
rep <- as.integer(args[6])
output_file <- args[7]
causal_snps_file <- args[8]

suppressMessages({
  library(data.table)
  library(dplyr)
})

message("=== Parameters ===")
message("Genotype pattern: ", geno_pattern)
message("GWAS pattern: ", gwas_pattern)
message("Test IDs: ", test_ids_file)
message("PLINK2 path: ", plink2_path)
message("Output dir: ", outdir)
message("Replicate: ", rep)
message("PRS output: ", output_file)
message("Causal SNPs file: ", causal_snps_file)

# Read test IDs and causal SNPs
test_ids <- fread(test_ids_file, header = FALSE)$V1
causal_snps <- fread(causal_snps_file, header = FALSE)$V1

# Initialize PRS vector (NA表示缺失)
prs_scores <- rep(NA_real_, length(test_ids))
names(prs_scores) <- test_ids

# Temporary directory for PLINK2 outputs
tmp_dir <- file.path(outdir, paste0("tmp_prs_rep", rep))
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over chromosomes
for(chr in 1:22){
  message("Processing chromosome ", chr)
  
  gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
  if(!file.exists(gwas_file)){
    message(" - GWAS file not found, skipping chr ", chr)
    next
  }
  
  # Filter GWAS to causal SNPs
  gwas_dt <- fread(gwas_file)
  gwas_dt <- gwas_dt %>% filter(ID %in% causal_snps)
  if(nrow(gwas_dt) == 0){
    message(" - No causal SNPs for chr ", chr)
    next
  }
  
  # Write filtered GWAS for PLINK2 scoring
  filtered_gwas <- file.path(tmp_dir, paste0("filtered_gwas_chr", chr, ".txt"))
  fwrite(gwas_dt, filtered_gwas, sep = "\t", quote = FALSE)
  
  # PLINK2 scoring
  plink_out <- file.path(tmp_dir, paste0("prs_chr", chr))
  system2(plink2_path, args = c(
    "--bfile", paste0(geno_pattern, chr, "_filtered"),
    "--score", filtered_gwas, "3 4 12 header cols=+scoresums",
    "--keep", test_ids_file,
    "--out", plink_out
  ))
  
  # Read .sscore
  sscore_file <- paste0(plink_out, ".sscore")
  if(!file.exists(sscore_file)){
    message(" - PRS score file missing for chr ", chr)
    next
  }
  
  sscore_dt <- fread(sscore_file)
  
  # Align individuals and accumulate PRS
  for(i in seq_len(nrow(sscore_dt))){
    iid <- sscore_dt$IID[i]
    idx <- match(iid, test_ids)
    if(!is.na(idx)){
      if(is.na(prs_scores[idx])){
        prs_scores[idx] <- sscore_dt$SCORE1_SUM[i]
      } else {
        prs_scores[idx] <- prs_scores[idx] + sscore_dt$SCORE1_SUM[i]
      }
    }
  }
  
  # Clean up temporary files for this chromosome
  file.remove(filtered_gwas)
  unlink(paste0(plink_out, c(".bed",".bim",".fam",".log")), recursive = TRUE)
}

# Ensure parent directory for output exists
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Write final PRS output
prs_df <- data.table(IID = test_ids, PRS_total = prs_scores)
fwrite(prs_df, output_file, sep = "\t")

# Clean up temporary directory
unlink(tmp_dir, recursive = TRUE)

message("PRS calculation complete. Output saved to: ", output_file)