#!/usr/bin/env Rscript

# ====== Config you can tweak (currently set for replicate 2) ======
rep <- 2
out_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/test"
ld_thresholds_str <- "0.5"  # add more like "0.1,0.2,0.5" if you clumped at those r^2
test_ids_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/test_ID_rep2.txt"
effects_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/test/effects_rep2.txt"
geno_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/ukb_chr"
plink_path <- "/data2/morgante_lab/yuxuans/software/plink2"

# Dominance term inclusion threshold (use p_d < p_threshold to include 'd_hat')
p_threshold <- 1e-3

suppressPackageStartupMessages({
  library(data.table)
  library(genio)
})

cat("\nStarting PRS calculation for replicate", rep, "\n")
cat("Output directory:", out_dir, "\n")
cat("LD thresholds:", ld_thresholds_str, "\n")
cat("Test ID file:", test_ids_file, "\n")
cat("Effect file:", effects_file, "\n")
cat("Geno pattern:", geno_pattern, "\n")
cat("PLINK path:", plink_path, "\n")
cat("Dominance p-threshold:", p_threshold, "\n")

# ------------------ helper: make clumped_snps files from PLINK2 .clumps ------------------
build_r2_clumped_snplist_if_missing <- function(out_dir, rep, r2_vec) {
  clump_files <- list.files(out_dir, pattern="\\.clumps$", full.names=TRUE)
  if (!length(clump_files)) {
    message("No .clumps files found in: ", out_dir, "  (skipping auto-build)")
    return(invisible(NULL))
  }
  one <- function(p) {
    dt <- tryCatch(fread(p), error=function(e) NULL)
    if (is.null(dt)) return(NULL)
    if ("#CHROM" %in% names(dt)) setnames(dt, "#CHROM", "CHR")
    if ("ID" %in% names(dt)) setnames(dt, "ID", "SNP")
    keep_cols <- intersect(c("CHR","SNP"), names(dt))
    if (!all(c("CHR","SNP") %in% keep_cols)) return(NULL)
    dt[, .(CHR=as.integer(CHR), SNP=as.character(SNP))]
  }
  combined <- rbindlist(lapply(clump_files, one), fill=TRUE)
  combined <- unique(combined[!is.na(CHR) & !is.na(SNP)])
  if (!nrow(combined)) {
    message("No usable rows in .clumps files.")
    return(invisible(NULL))
  }
  for (r2 in r2_vec) {
    out_path <- file.path(out_dir, sprintf("clumped_snps_rep%d_r2_%.1f.txt", rep, as.numeric(r2)))
    if (!file.exists(out_path)) {
      fwrite(combined, out_path, sep="\t")
      message("Wrote index SNP list for r2=", r2, " -> ", out_path)
    } else {
      message("Found existing: ", out_path)
    }
  }
}

# ------------------ temp dir ------------------
temp_dir <- file.path(out_dir, sprintf("temp_rep%d", rep))
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created temp directory:", temp_dir, "\n")
}

# ------------------ thresholds & inputs ------------------
ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds_str, ",")))
cat("LD thresholds parsed:", paste(ld_thresholds, collapse=", "), "\n")

test_ids <- fread(test_ids_file, header = FALSE)[[1]]
cat("✅ Loaded", length(test_ids), "test individuals\n")

eff <- fread(effects_file)
if ("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")
if ("#CHROM" %in% names(eff)) setnames(eff, "#CHROM", "CHR")
missing_cols <- setdiff(c("CHR","SNP","ref_allele","alt_allele","a_hat","d_hat","p_a","p_d"), names(eff))
if (length(missing_cols)) stop("Effects table missing columns: ", paste(missing_cols, collapse=", "))

# Build clumped SNP lists from .clumps (once per r2)
build_r2_clumped_snplist_if_missing(out_dir, rep, ld_thresholds)

# ------------------ helpers ------------------
safe_add <- function(x, y) {
  out <- x
  both <- !is.na(x) & !is.na(y)
  only_y <-  is.na(x) & !is.na(y)
  out[both] <- x[both] + y[both]
  out[only_y] <- y[only_y]
  out
}

# Fast per-SNP (column-wise) mode imputation for genotypes (0/1/2)
impute_mode_matrix <- function(X) {
  # counts for each genotype; NAs ignored
  c0 <- colSums(X == 0, na.rm = TRUE)
  c1 <- colSums(X == 1, na.rm = TRUE)
  c2 <- colSums(X == 2, na.rm = TRUE)
  # mode code per column; ties resolved by preferring 0, then 1, then 2
  mcode <- apply(cbind(c0, c1, c2), 1, function(cnt) which.max(cnt) - 1)
  na_idx <- which(is.na(X), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    X[na_idx] <- mcode[na_idx[, "col"]]
  }
  X
}

# ------------------ PRS storage ------------------
ncol_prs <- length(ld_thresholds)
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = ncol_prs,
                     dimnames = list(test_ids, paste0("r2_", format(ld_thresholds, nsmall = 1))))

# ------------------ main loop over r2 ------------------
for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  cat("\n=== LD threshold r² =", r2_str, "===\n")
  
  prs_vector <- numeric(length(test_ids))
  prs_vector[] <- 0
  
  clump_file <- file.path(out_dir, sprintf("clumped_snps_rep%d_r2_%.1f.txt", rep, r2))
  cat("Clump SNP list:", clump_file, "\n")
  if (!file.exists(clump_file)) {
    cat("  ⚠️  Clump SNP list not found, skipping r2 =", r2_str, "\n")
    next
  }
  clumped <- fread(clump_file)
  if ("#CHROM" %in% names(clumped)) setnames(clumped, "#CHROM", "CHR")
  if ("ID" %in% names(clumped)) setnames(clumped, "ID", "SNP")
  if (!all(c("CHR","SNP") %in% names(clumped))) stop("Clump SNP list must contain CHR and SNP columns.")
  
  # Join effects for these clumped SNPs
  clumped_effects <- merge(
    unique(clumped[, .(SNP, CHR = as.integer(CHR))]),
    eff[, .(SNP, CHR, ref_allele, alt_allele, a_hat, d_hat, p_a, p_d)],
    by = c("SNP", "CHR")
  )
  if (nrow(clumped_effects) == 0) {
    cat("  No overlapping SNPs between clumped list and effects for r2 =", r2_str, "\n")
    prs_matrix[, r2_idx] <- prs_vector
    next
  }
  
  # -------------- per chromosome --------------
  for (chr in 1:22) {
    cat("  Processing chromosome", chr, "...\n")
    chr_snps <- clumped_effects[CHR == chr,
                                .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_a, p_d)]
    if (nrow(chr_snps) == 0) {
      cat("    No SNPs for chromosome", chr, "\n")
      next
    }
    
    # Extract SNPs for test set
    snp_list_file <- file.path(temp_dir, sprintf("chr%d_snps.txt", chr))
    fwrite(chr_snps[, .(SNP)], snp_list_file, col.names = FALSE)
    
    plink_file <- paste0(geno_pattern, chr, "_filtered")
    temp_prefix <- file.path(temp_dir, sprintf("chr%d_r2_%s", chr, r2_str))
    
    plink_cmd <- sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s --memory 2000 --threads 10",
      shQuote(plink_path),
      shQuote(plink_file),
      shQuote(test_ids_file),
      shQuote(snp_list_file),
      shQuote(temp_prefix)
    )
    system(plink_cmd)
    
    bed_file <- paste0(temp_prefix, ".bed")
    bim_file <- paste0(temp_prefix, ".bim")
    fam_file <- paste0(temp_prefix, ".fam")
    if (!file.exists(bed_file) || !file.exists(bim_file) || !file.exists(fam_file)) {
      cat("    No extracted genotype data for chr", chr, "\n")
      next
    }
    
    plink_data <- genio::read_plink(temp_prefix, verbose = FALSE)
    fam_ids <- plink_data$fam$id
    snp_ids <- plink_data$bim$id
    X <- t(plink_data$X)  # Individuals x SNPs (0/1/2/NA)
    
    # Reorder individuals to match test_ids
    if (!all(fam_ids == test_ids)) {
      idx_map <- match(test_ids, fam_ids)
      if (any(is.na(idx_map))) stop("Mismatch between test_ids and FAM IDs after --keep.")
      X <- X[idx_map, , drop = FALSE]
      fam_ids <- fam_ids[idx_map]
    }
    
    # Align effect rows to SNP order in X
    snp_effects <- chr_snps[match(snp_ids, chr_snps$SNP), ]
    if (any(is.na(snp_effects$SNP))) {
      keep <- !is.na(snp_effects$SNP)
      X <- X[, keep, drop = FALSE]
      snp_ids <- snp_ids[keep]
      snp_effects <- snp_effects[keep]
    }
    if (!length(snp_ids)) {
      cat("    After alignment, 0 SNPs remain for chr", chr, "\n")
      next
    }
    
    # -------- Mode imputation (per SNP/column) --------
    # Replace NA genotypes with modal value among {0,1,2} for that SNP
    X <- impute_mode_matrix(X)
    
    # Accumulate PRS
    prs_contrib <- numeric(nrow(X))
    for (j in seq_along(snp_ids)) {
      g <- X[, j]  # 0/1/2 (no NA after imputation)
      a <- snp_effects$a_hat[j]
      d <- snp_effects$d_hat[j]
      p_d_val <- snp_effects$p_d[j]
      
      contrib <- numeric(length(g))
      if (!is.na(p_d_val) && p_d_val < p_threshold) {
        # include dominance
        contrib[g == 0] <- -a
        contrib[g == 1] <-  d
        contrib[g == 2] <-  a
      } else {
        # additive only
        contrib[g == 0] <- -a
        contrib[g == 1] <-  0
        contrib[g == 2] <-  a
      }
      prs_contrib <- safe_add(prs_contrib, contrib)
    }
    
    prs_vector <- safe_add(prs_vector, prs_contrib)
  } # end chr loop
  
  prs_matrix[, r2_idx] <- prs_vector
} # end r2 loop

# Cleanup temp
unlink(temp_dir, recursive = TRUE)

# Save PRS
prs_df <- data.table(IID = test_ids, prs_matrix)
output_file <- file.path(out_dir, sprintf("test_prs_scan_impute_rep%d.txt", rep))
fwrite(prs_df, output_file, sep = "\t")
cat("✅ Saved PRS results to:", output_file, "\n")

