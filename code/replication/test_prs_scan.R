#!/usr/bin/env Rscript

# === Manual PRS Calculation for Scenario 4, rep=1, p=10E-5, LD=0.1 ===

# Parameters
rep <- 1
p_threshold <- "10E-5"  # Used as string in file paths and for logging
ld_thresholds_str <- "0.1"
out_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/10E-5"
temp_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/test_prs_gwas"
test_ids_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/test_ID_rep1.txt"
geno_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/ukb_chr"
gwas_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/gwas/gwas_rep1_chr"
clump_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/clumping/10E-5/gwas_clump_rep1_chr"
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
output_file <- file.path(temp_dir, "temp_prs_gwas_rep1.txt")

# Derived parameters
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))

# Create temp directory
temp_dir <- file.path(temp_dir, sprintf("temp_prs_rep%d", rep))
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
  cat("Created temporary directory:", temp_dir, "\n")
}

# Load required package
library(data.table)

test_ids <- fread(test_ids_file, header = FALSE)$V1
cat("✅ Loaded", length(test_ids), "test individuals\n")

# Initialize PRS matrix
ncol_prs <- length(ld_thresholds)
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = ncol_prs,
                     dimnames = list(test_ids, paste0("r2_", format(ld_thresholds, nsmall = 1))))

# Process chromosomes 1-22
for (chr in 1:22) {
  cat("🧬 Chromosome", chr, "processing... \n")
  
  # File paths
  gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
  geno_prefix <- paste0(geno_pattern, chr, "_filtered")
  temp_prefix <- file.path(temp_dir, sprintf("chr%d", chr))
  
  # Read GWAS results
  gwas <- fread(gwas_file)
  
  # Process each LD threshold
  for (r2 in ld_thresholds) {
    r2_str <- format(r2, nsmall = 1)
    clump_file <- paste0(clump_pattern, chr, "_r2_", r2_str, ".snplist")
    score_file <- paste0(temp_prefix, "_r2_", r2_str, ".score")
    extract_prefix <- paste0(temp_prefix, "_r2_", r2_str)
    
    
    # Read clumped SNPs (single column)
    clumped_snps <- fread(clump_file, header = FALSE)
    
    # Filter GWAS to clumped SNPs
    weights <- gwas[ID %in% clumped_snps$V1, .(ID, REF, ALT, A1, BETA)]
    
    # Set A2 = REF
    weights[, A2 := REF]
    weights[, BETA := as.numeric(BETA)]
    
    # Prepare score file: SNP, A1 (effect allele), WEIGHT
    setnames(weights, c("ID", "BETA"), c("SNP", "WEIGHT"))
    weights <- weights[, .(SNP, A1, WEIGHT)]
    
    # Write score file
    fwrite(weights, score_file, sep = "\t", quote = FALSE, col.names = FALSE)
    
    # Extract genotypes for test individuals
    plink_extract_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s --silent --memory 2000 --threads 40",
      plink2_path, geno_prefix, test_ids_file, clump_file, extract_prefix
    )
    cat("  Running PLINK2 extract...\n")
    system(plink_extract_cmd)
    
    # Compute PRS
    plink_score_cmd <- sprintf(
      "%s --bfile %s --score %s cols=+scoresums --out %s --silent",
      plink2_path, extract_prefix, score_file, extract_prefix
    )
    cat("  Running PLINK2 score...\n")
    system(plink_score_cmd)
    
    # Read PRS chunk
    prs_chunk_file <- paste0(extract_prefix, ".sscore")
    
    prs_chunk <- fread(prs_chunk_file)
    
    # Match sample IDs
    idx <- match(prs_chunk$IID, test_ids)
    valid <- !is.na(idx)
    
    # Add to PRS matrix
    col_idx <- which(ld_thresholds == r2)
    prs_matrix[idx[valid], col_idx] <- prs_matrix[idx[valid], col_idx] + prs_chunk$SCORE1_SUM[valid]
    
    # Optional cleanup
    unlink(paste0(extract_prefix, c(".bed", ".bim", ".fam", ".log", ".nosex", ".score", ".sscore")), force = TRUE)
  }
  
  cat("✅ Done (chromosome", chr, ")\n")
}

# Save final PRS
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

prs_df <- data.table(IID = test_ids, prs_matrix)

fwrite(prs_df, output_file, sep = "\t")
cat("\n🎉 SUCCESS: Full PRS computed and saved to:\n")
cat("   ", output_file, "\n")
cat("   Total individuals:", nrow(prs_df), "\n")
cat("   PRS range:", range(prs_df$r2_0.1), "\n")

# Optional: Show first few lines
cat("\n📊 First 10 PRS values:\n")
print(prs_df[1:10])

# Final cleanup
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  cat("🧹 Temporary directory removed.\n")
}

cat("✅ Full 22-chromosome PRS calculation completed successfully!\n")

pheno <- fread("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/simulated_pheno_rep1.txt")
merged = merge(prs_df, pheno, by = "IID", all.x = T)
cor(merged$r2_0.1, merged$pheno)
cor(merged$r2_0.1, merged$pheno)^2
