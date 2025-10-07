#!/usr/bin/env Rscript

# ---------------- Parse args ----------------
args <- commandArgs(trailingOnly = TRUE)

geno_pattern    <- args[1]
gwas_pattern    <- args[2]
clump_pattern   <- args[3]
test_ids_file   <- args[4]
plink2_path     <- args[5]
ld_thresholds   <- as.numeric(unlist(strsplit(args[6], ",")))
out_dir         <- args[7]
rep             <- args[8]
p_threshold     <- args[9]
output_file     <- args[10]

cat("\n[PRS GWAS] Replicate", rep, " P-threshold", p_threshold, "\n")

# ---------------- Setup ----------------
suppressPackageStartupMessages({
  library(data.table)
})

if (!file.exists(test_ids_file)) stop("❌ Missing test_ids file: ", test_ids_file)
test_ids <- fread(test_ids_file, header = FALSE)[[1]]

prs_matrix <- matrix(NA_real_, nrow = length(test_ids), ncol = length(ld_thresholds),
                     dimnames = list(test_ids, paste0("r2_", format(ld_thresholds, nsmall = 1))))

temp_dir <- file.path(out_dir, paste0("temp_prs_rep", rep))
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Loop chromosomes ----------------
for (chr in 1:22) {
  cat("🧬 Chr", chr, "...\n")
  gwas_file   <- paste0(gwas_pattern, chr, ".pheno.glm.linear")
  geno_prefix <- paste0(geno_pattern, chr, "_filtered")
  temp_prefix <- file.path(temp_dir, paste0("chr", chr))
  
  if (!file.exists(gwas_file) || !file.exists(paste0(geno_prefix, ".bim"))) {
    cat("  ↳ GWAS/geno missing, skip\n")
    next
  }
  
  gwas <- fread(gwas_file)
  if (!("ID" %in% names(gwas)) && "SNP" %in% names(gwas)) setnames(gwas, "SNP", "ID")
  gwas_add <- gwas[TEST == "ADD"]   # ✅ 只用 ADD 行
  
  for (r2 in ld_thresholds) {
    r2_str <- format(r2, nsmall = 1)
    col_idx <- which(ld_thresholds == r2)
    
    clump_file   <- paste0(clump_pattern, chr, "_r2_", r2_str, ".snplist")
    extract_pref <- paste0(temp_prefix, "_r2_", r2_str)
    score_file   <- paste0(extract_pref, ".score")
    sscore_file  <- paste0(extract_pref, ".sscore")
    
    if (!file.exists(clump_file)) { next }
    clumped_snps <- fread(clump_file, header = FALSE)
    if (nrow(clumped_snps) == 0L) { next }
    setnames(clumped_snps, "V1", "ID")
    
    # ✅ 只筛选 ADD + 去重
    weights <- unique(
      gwas_add[ID %in% clumped_snps$ID, .(SNP = ID, A1, WEIGHT = as.numeric(BETA))],
      by = "SNP"
    )
    weights <- weights[!is.na(WEIGHT) & !is.na(A1) & !is.na(SNP)]
    if (nrow(weights) == 0L) { next }
    
    fwrite(weights, score_file, sep = "\t", quote = FALSE, col.names = FALSE)
    
    # --- Extract test genotypes ---
    system(sprintf(
      '%s --bfile %s --keep %s --extract %s --make-bed --out %s --silent',
      shQuote(plink2_path), shQuote(geno_prefix), shQuote(test_ids_file),
      shQuote(clump_file), shQuote(extract_pref)
    ))
    if (!file.exists(paste0(extract_pref, ".bim"))) { next }
    
    # --- Compute PRS ---
    system(sprintf(
      '%s --bfile %s --score %s cols=+scoresums list-variants --out %s --silent',
      shQuote(plink2_path), shQuote(extract_pref), shQuote(score_file), shQuote(extract_pref)
    ))
    if (!file.exists(sscore_file)) { next }
    
    prs_chunk <- fread(sscore_file)
    if (!"IID" %in% names(prs_chunk)) { next }
    
    score_col <- if ("SCORE1_SUM" %in% names(prs_chunk)) "SCORE1_SUM"
    else grep("^SCORE.*_SUM$", names(prs_chunk), value = TRUE)[1]
    if (is.na(score_col)) { next }
    
    idx   <- match(prs_chunk$IID, test_ids)
    valid <- !is.na(idx)
    new_scores <- suppressWarnings(as.numeric(prs_chunk[[score_col]]))
    
    current <- prs_matrix[idx[valid], col_idx]
    addend  <- new_scores[valid]
    updated <- ifelse(is.na(addend), current,
                      ifelse(is.na(current), addend, current + addend))
    prs_matrix[idx[valid], col_idx] <- updated
    
    # cleanup
    unlink(paste0(extract_pref, c(".bed", ".bim", ".fam", ".log", ".nosex",
                                  ".sscore", ".sscore.vars")), force = TRUE)
    unlink(score_file, force = TRUE)
  }
}

# ---------------- Save ----------------
prs_df <- data.table(IID = test_ids, prs_matrix)
fwrite(prs_df, output_file, sep = "\t")
cat("✅ Saved:", output_file, " with", nrow(prs_df), "samples\n")
