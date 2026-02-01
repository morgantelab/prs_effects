#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(genio)
})

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# Args (STRICT ORDER)
# -----------------------------
effects_file   <- args[1]   # ad_clumped.txt (CHR POS SNP REF ALT a_hat d_hat p_a p_d)
bfile_merged   <- args[2]   # merged bfile prefix
test_ids_file  <- args[3]   # keep file (FID IID)
pheno_file     <- args[4]   # PHENO_CSV
trait_name     <- args[5]   # e.g., MCV
out_prs        <- args[6]
out_corr       <- args[7]
plink2_path    <- args[8]
dom_p_thresh   <- as.numeric(args[9])
threads_plink  <- as.integer(args[10])

# -----------------------------
# Helpers
# -----------------------------
impute_mode_matrix <- function(X) {
  # X: numeric/integer matrix with NA, individuals x SNPs
  if (!is.matrix(X) || ncol(X) == 0) return(X)
  for (j in seq_len(ncol(X))) {
    g <- X[, j]
    if (anyNA(g)) {
      gg <- g[!is.na(g)]
      if (length(gg) == 0L) {
        X[, j] <- 0L
      } else {
        tab <- tabulate(gg + 1L, nbins = 3L)  # assumes values in 0/1/2
        mode_val <- which.max(tab) - 1L
        g[is.na(g)] <- mode_val
        X[, j] <- g
      }
    }
  }
  X
}

# -----------------------------
# Load inputs
# -----------------------------
effects <- fread(effects_file)

req_cols <- c("CHR","SNP","REF","ALT","a_hat","d_hat","p_d")
if (!all(req_cols %in% names(effects))) {
  stop("effects table missing required columns. Need: ", paste(req_cols, collapse = ", "))
}
effects[, CHR := as.integer(CHR)]

test_ids_dt  <- fread(test_ids_file, header = FALSE)
setnames(test_ids_dt, c("FID","IID"))
test_ids_vec <- test_ids_dt[["IID"]]

ph <- fread(pheno_file)
if (!(trait_name %in% names(ph))) {
  stop("Trait column not found in phenotype file: ", trait_name)
}
pheno <- ph[, .(IID, pheno = get(trait_name))]

# -----------------------------
# PRS accumulation (Method2 logic, using merged bfile)
# -----------------------------
tmp_dir <- file.path(dirname(out_prs), "tmp_prs_ad_merged")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

prs_vec <- numeric(length(test_ids_vec))
names(prs_vec) <- test_ids_vec

for (chr in sort(unique(effects$CHR))) {
  se_chr <- effects[CHR == chr, .(SNP, REF, ALT, a_hat, d_hat, p_d)]
  if (!nrow(se_chr)) next

  # write snplist for this chr
  snplist_chr <- file.path(tmp_dir, sprintf("prs_chr%d.snps", chr))
  fwrite(se_chr[, .(SNP)], snplist_chr, col.names = FALSE)

  # extract genotypes for these SNPs from merged bfile
  out_prefix <- file.path(tmp_dir, sprintf("extract_chr%d", chr))
  cmd <- sprintf(
    "%s --bfile %s --keep %s --extract %s --make-bed --out %s --threads %d",
    shQuote(plink2_path),
    shQuote(bfile_merged),
    shQuote(test_ids_file),
    shQuote(snplist_chr),
    shQuote(out_prefix),
    threads_plink
  )

  ret <- system(cmd)
  if (ret != 0L || !file.exists(paste0(out_prefix, ".bed"))) next

  gd <- genio::read_plink(out_prefix, verbose = FALSE)

  # genio returns X as SNP x individual -> transpose to individual x SNP
  X <- t(gd$X)

  # Align rows to test IDs
  fam_ids <- gd$fam$id
  idx <- match(test_ids_vec, fam_ids)
  X <- X[idx, , drop = FALSE]

  # Align columns to effects order
  snp_ids <- gd$bim$id
  se <- se_chr[match(snp_ids, se_chr$SNP)]
  keep <- !is.na(se$SNP)
  X <- X[, keep, drop = FALSE]
  se <- se[keep]

  if (ncol(X) == 0) next

  # Impute missing genotypes (mode)
  X <- impute_mode_matrix(X)

  # Add contributions SNP-by-SNP (vectorized per SNP)
  for (j in seq_len(ncol(X))) {
    g  <- X[, j]
    a  <- se$a_hat[j]
    d  <- se$d_hat[j]
    pd <- se$p_d[j]

    # additive: 0->-a, 1->0, 2->+a
    contrib <- (g - 1) * a

    # dominance only if p_d < threshold, and only for het (g==1)
    if (!is.na(pd) && pd < dom_p_thresh) {
      contrib <- contrib + (g == 1L) * d
    }

    prs_vec <- prs_vec + contrib
  }
}

# -----------------------------
# Save PRS + correlation
# -----------------------------
dir.create(dirname(out_prs), showWarnings = FALSE, recursive = TRUE)
fwrite(data.table(IID = test_ids_vec, PRS = prs_vec), out_prs, sep = "\t")

prs_dt <- data.table(IID = test_ids_vec, PRS = prs_vec)
prs_dt[, IID := as.character(IID)]

pheno[, IID := as.character(IID)]
res <- merge(prs_dt, pheno, by = "IID")
rval <- suppressWarnings(cor(res$pheno, res$PRS, use = "complete.obs"))

dir.create(dirname(out_corr), showWarnings = FALSE, recursive = TRUE)
fwrite(data.table(corr_ad = rval, N = nrow(res)), out_corr, sep = "\t")

cat(sprintf("[TEST merged method2 a+d] trait=%s N=%d cor=%.6f\n", trait_name, nrow(res), rval))
cat("Saved PRS -> ", out_prs, "\n", sep = "")
