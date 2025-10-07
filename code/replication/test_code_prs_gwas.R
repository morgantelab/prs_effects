# ===================================
# FULL PRS CALCULATION SCRIPT
# For: rep=1, p_threshold=10E-3, r2=0.1
# Processes all 22 chromosomes
# ===================================

# === 1. Hardcoded Paths from Config ===
pheno_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3/pheno"
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3"
code_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/code/replication"

plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

# === 2. Run Parameters ===
rep <- "1"
p_threshold <- "10E-3"
r2 <- 0.1
r2_str <- format(r2, nsmall = 1)  # "0.1"

# === 3. File Patterns ===
geno_pattern <- file.path(pheno_dir, "ukb_chr")  # → ukb_chr1_filtered.bed
gwas_pattern <- file.path(results_dir, "gwas", sprintf("gwas_rep%s_chr", rep))  # gwas_rep1_chr1.pheno.glm.linear
clump_pattern <- file.path(results_dir, "clumping", p_threshold, 
                           sprintf("gwas_clump_rep%s_chr", rep))  # gwas_clump_rep1_chr1_r2_0.1.snplist
test_ids_file <- file.path(pheno_dir, sprintf("test_ID_rep%s.txt", rep))

# Output
out_dir <- file.path(results_dir, "test", p_threshold)
output_file <- file.path(out_dir, sprintf("test_prs_gwas_rep%s_r2%s.txt", rep, r2_str))
temp_dir <- file.path(out_dir, paste0("temp_prs_rep", rep))

# Create directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# === 4. Load Required Package ===
library(data.table)

# === 5. Load Test IDs ===
if (!file.exists(test_ids_file)) {
  stop("❌ Test IDs file not found: ", test_ids_file)
}
test_ids <-fread(test_ids_file, header = FALSE)$V1

# === 6. Initialize PRS Matrix (one column for r² = 0.1) ===
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = 1,
                     dimnames = list(test_ids, paste0("r2_", r2_str)))

# === 7. Loop Over Chromosomes 1 to 22 ===
for (chr in 1:22) {
  cat("🧬 Chromosome", chr, "processing... ", "\n")
  
  # File paths
  clump_file <- paste0(clump_pattern, chr, "_r2_", r2_str, ".snplist")
  gwas_file <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
  geno_prefix <- paste0(geno_pattern, chr, "_filtered")
  temp_prefix <- file.path(temp_dir, sprintf("chr%d_r2%s", chr, r2_str))
  score_file <- paste0(temp_prefix, ".score")
  profile_file <- paste0(temp_prefix, ".sscore")
  
  # --- Step 1: Check Input Files ---
  if (!file.exists(clump_file)) {
    cat("❌ Clump file missing:", clump_file, "\n")
    next
  }
  if (!file.exists(gwas_file)) {
    cat("❌ GWAS file missing:", gwas_file, "\n")
    next
  }
  if (!file.exists(paste0(geno_prefix, ".bim"))) {
    cat("❌ Genotype .bim file missing for chr", chr, "\n")
    next
  }
  
  # --- Step 2: Read Clumped SNPs ---
  clumped_snps <- fread(clump_file, header = FALSE)[[1]]
  
  # --- Step 3: Read GWAS Results ---
  gwas <-fread(gwas_file)

  # Filter to clumped SNPs
  weights <- gwas[ID %in% clumped_snps, .(ID, REF, ALT, A1, BETA)]
  
  # With --glm omit-ref: A1 = ALT (effect), A2 = REF (baseline)
  weights[, A2 := REF]
  weights[, BETA := as.numeric(BETA)]
  
  # Prepare score file for PLINK2
  setnames(weights, c("ID", "BETA"), c("SNP", "WEIGHT"))
  setcolorder(weights, c("SNP", "A1", "A2", "WEIGHT"))
  weights <- weights %>% select(SNP, A1, WEIGHT)
  
  fwrite(weights, score_file, sep = "\t", quote = FALSE, col.names = FALSE)
  
  # --- Step 4: Extract Genotypes with PLINK2 ---
  cat("Extracting SNPs for chromosome", chr,  "\n")
  plink_extract_cmd <- sprintf(
    "%s --bfile %s --keep %s --extract %s --make-bed --out %s --silent",
    plink2_path, geno_prefix, test_ids_file, clump_file, temp_prefix
  )
  system(plink_extract_cmd)
  
  
  # --- Step 5: Compute PRS with PLINK2 ---
  cat("Computing PRS for chromosome", chr,  "\n")
  
  plink_score_cmd <- sprintf(
    "%s --bfile %s --score %s cols=+scoresums --out %s --silent",
    plink2_path, temp_prefix, score_file, temp_prefix
  )
  system(plink_score_cmd)
  
  # --- Step 6: Read and Accumulate PRS ---
  prs_chunk <- fread(profile_file)
  
  # Match sample IDs
  idx <- match(prs_chunk$IID, test_ids)
  valid <- !is.na(idx)
  
  # Add to PRS matrix
  prs_matrix[idx[valid], 1] <- prs_matrix[idx[valid], 1] + prs_chunk$SCORE1_SUM[valid]
  
  # Cleanup
  # unlink(c(
  #   paste0(temp_prefix, ".bed"),
  #   paste0(temp_prefix, ".bim"),
  #   paste0(temp_prefix, ".fam"),
  #   paste0(temp_prefix, ".log"),
  #   paste0(temp_prefix, ".nosex"),
  #   paste0(temp_prefix, ".score"),
  #   paste0(temp_prefix, ".profile"),
  #   paste0(temp_prefix, ".missing")
  # ), force = TRUE)
  # 
  cat("Done (", nrow(weights), "SNPs)\n")
}

# === 8. Save Final PRS ===
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

prs_df <- data.table(IID = test_ids, prs_matrix)
fwrite(prs_df, output_file, sep = "\t")
cat("\n🎉 SUCCESS: Full PRS computed and saved to:\n")
cat("   ", output_file, "\n")
cat("   Total individuals:", nrow(prs_df), "\n")
cat("   PRS range:", range(prs_df[[2]]), "\n")

# Optional: Show first few lines
cat("\n📊 First 10 PRS values:\n")
print(prs_df[1:10])

# === 9. Final Cleanup (Uncomment when confident)
if (dir.exists(temp_dir)) {
  unlink(temp_dir, recursive = TRUE)
  cat("🧹 Temporary directory removed.\n")
}

cat("✅ Full 22-chromosome PRS calculation completed successfully!\n")