#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
  library(genio)
})

# -------------------- Args --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop("Usage: estimate_effects_and_prs_ad.R ",
       "PHENO_DIR PLINK2 FILTERED_PREFIX TRAIN_IDS TEST_IDS TRAIN_PHENO TEST_PHENO ",
       "DOM_PTHRESH CHUNK_SIZE THREADS NCHR")
}

pheno_dir       <- args[1]   # e.g. /.../OUT_DIR/<scen>/rep<rep>
plink2_path     <- args[2]   # e.g. /path/to/plink2
filtered_prefix <- args[3]   # e.g. /.../OUT_DIR/filtered/ukb_chr
train_ids_file  <- args[4]
test_ids_file   <- args[5]
train_pheno_file<- args[6]
test_pheno_file <- args[7]
dom_p_thresh    <- as.numeric(args[8])
chunk_size      <- as.integer(args[9])
threads_plink   <- as.integer(args[10])
NCHR            <- as.integer(args[11])

stopifnot(file.exists(plink2_path))

# Files generated previously and expected here
snplist_file <- file.path(pheno_dir, "clump", "clumped_snps.snplist")
if (!file.exists(snplist_file)) {
  stop("Missing clumped SNP list: ", snplist_file)
}

# Load IDs and phenotypes
train_ids <- fread(train_ids_file, header = FALSE, col.names = c("FID","IID"))
test_ids_df <- fread(test_ids_file, header = FALSE, col.names = c("FID","IID"))

training_pheno <- fread(
  train_pheno_file,
  colClasses = c(FID="character", IID="character", pheno="numeric")
)
training_pheno[, IID := as.character(IID)]

test_pheno <- fread(test_pheno_file)
test_pheno[, IID := as.character(IID)]

# -------------------- Helpers --------------------
calculate_welch <- function(data, snp_cols) {
  n_snps <- length(snp_cols)
  results <- data.frame(
    SNP   = snp_cols,
    a_hat = rep(NA_real_, n_snps),
    d_hat = rep(NA_real_, n_snps),
    p_a   = rep(NA_real_, n_snps),
    p_d   = rep(NA_real_, n_snps)
  )

  pheno <- data$pheno

  for (j in seq_len(n_snps)) {
    snp <- snp_cols[j]
    X   <- data[[snp]]

    AA <- X == 2L
    Aa <- X == 1L
    aa <- X == 0L

    n_AA <- sum(AA, na.rm = TRUE)
    n_Aa <- sum(Aa, na.rm = TRUE)
    n_aa <- sum(aa, na.rm = TRUE)

    if (n_AA < 2 || n_Aa < 2 || n_aa < 2) next

    y_AA <- pheno[AA]
    y_Aa <- pheno[Aa]
    y_aa <- pheno[aa]

    mean_AA <- mean(y_AA, na.rm = TRUE)
    mean_Aa <- mean(y_Aa, na.rm = TRUE)
    mean_aa <- mean(y_aa, na.rm = TRUE)
    var_AA  <- var(y_AA,  na.rm = TRUE)
    var_Aa  <- var(y_Aa,  na.rm = TRUE)
    var_aa  <- var(y_aa,  na.rm = TRUE)

    # Additive (Welch)
    a_hat <- (mean_AA - mean_aa) / 2
    var_a <- 0.25 * (var_AA / n_AA + var_aa / n_aa)
    if (!is.finite(var_a) || var_a <= 0) next

    df_a <- (var_AA/n_AA + var_aa/n_aa)^2 /
      ((var_AA/n_AA)^2/(n_AA-1) + (var_aa/n_aa)^2/(n_aa-1))
    t_a <- a_hat / sqrt(var_a)
    p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)

    # Dominance (Welch)
    d_hat <- mean_Aa - (mean_AA + mean_aa)/2
    var_d <- var_Aa/n_Aa + 0.25*(var_AA/n_AA + var_aa/n_aa)
    if (!is.finite(var_d) || var_d <= 0) next

    df_d <- (var_d)^2 / (
      (var_Aa/n_Aa)^2/(n_Aa-1) +
        (var_AA/(4*n_AA))^2/(n_AA-1) +
        (var_aa/(4*n_aa))^2/(n_aa-1)
    )
    t_d <- d_hat / sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)

    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  results
}

impute_mode_matrix <- function(X) {
  nc <- ncol(X)
  for (j in seq_len(nc)) {
    if (anyNA(X[, j])) {
      g <- X[, j]
      tab <- tabulate(g[!is.na(g)] + 1L, nbins = 3L)
      g[is.na(g)] <- which.max(tab) - 1L
      X[, j] <- g
    }
  }
  X
}

# -------------------- 1) Estimate a,d on TRAINING (Welch) --------------------
all_results <- list()

for (chr in seq_len(NCHR)) {
  message("Estimating effects on chromosome ", chr)

  full_plink <- sprintf("%s%d_filtered", filtered_prefix, chr)
  if (!file.exists(paste0(full_plink, ".bim"))) {
    message("  - No BIM for chr ", chr, " -> skip")
    next
  }

  tmp_prefix <- file.path(tempdir(), sprintf("chr%d_train_ct", chr))

  plink_cmd <- sprintf(
    '%s --bfile %s --keep %s --keep-allele-order --extract %s --make-bed --out %s --threads %d',
    shQuote(plink2_path),
    shQuote(full_plink),
    shQuote(train_ids_file),
    shQuote(snplist_file),
    shQuote(tmp_prefix),
    threads_plink
  )
  message("  Running: ", plink_cmd)
  ret <- system(plink_cmd)
  if (ret != 0L) {
    message("  - plink2 extract failed on chr ", chr, ", skipping")
    next
  }

  tmp_bed <- paste0(tmp_prefix, ".bed")
  tmp_bim <- paste0(tmp_prefix, ".bim")
  tmp_fam <- paste0(tmp_prefix, ".fam")
  if (!file.exists(tmp_bed)) {
    message("  - Empty extract on chr ", chr, " -> skip")
    next
  }

  fam <- fread(tmp_fam, col.names = c("FID","IID","pat","mat","sex","pheno"))
  bim <- fread(tmp_bim, col.names = c("chr","SNP","pos_cM","pos","alt","ref"))

  pheno_vec <- training_pheno$pheno[match(fam$IID, training_pheno$IID)]
  G <- BEDMatrix(tmp_bed)
  n_snps_total <- ncol(G)
  if (is.null(n_snps_total) || n_snps_total == 0) {
    unlink(paste0(tmp_prefix, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
    next
  }

  chr_results <- vector("list", ceiling(n_snps_total / chunk_size))
  idx <- 1L
  for (start_col in seq(1L, n_snps_total, by = chunk_size)) {
    end_col <- min(start_col + chunk_size - 1L, n_snps_total)
    snp_idx <- start_col:end_col

    snp_ids_block <- bim$SNP[snp_idx]
    snp_block <- G[, snp_idx, drop = FALSE]
    colnames(snp_block) <- snp_ids_block

    df_block <- data.frame(
      pheno = pheno_vec,
      as.data.frame(snp_block),
      check.names = FALSE
    )
    res_block <- calculate_welch(df_block, snp_ids_block)
    chr_results[[idx]] <- res_block

    message(sprintf("    chr %d: %d-%d / %d", chr, start_col, end_col, n_snps_total))
    idx <- idx + 1L
  }

  results_dt <- rbindlist(chr_results, use.names = TRUE, fill = TRUE)
  results_final <- merge(
    results_dt,
    bim[, .(SNP, pos, ref, alt)],
    by = "SNP",
    all.x = TRUE
  )[, .(
    chromosome = chr,
    position   = pos,
    SNP,
    ref_allele = ref,
    alt_allele = alt,
    a_hat, d_hat, p_a, p_d
  )]

  all_results[[length(all_results) + 1L]] <- results_final
  unlink(paste0(tmp_prefix, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
}

effects_out <- file.path(pheno_dir, "effects", "effects_training_method2.txt")
dir.create(dirname(effects_out), showWarnings = FALSE, recursive = TRUE)
fwrite(rbindlist(all_results, use.names = TRUE, fill = TRUE), effects_out, sep = "\t")
cat("✅ Wrote effects → ", effects_out, "\n", sep = "")

# -------------------- 2) PRS a+d on TEST --------------------
test_ids_vec <- fread(test_ids_file, header = FALSE)[[2]]

tmp_dir <- file.path(pheno_dir, "tmp_prs_ad")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

effects <- fread(effects_out)
if (!all(c("SNP","a_hat","d_hat","p_d","chromosome") %in% names(effects))) {
  stop("effects table missing required columns")
}
effects[, chromosome := as.integer(chromosome)]

prs_vec <- numeric(length(test_ids_vec)); names(prs_vec) <- test_ids_vec

for (chr in seq_len(NCHR)) {
  se_chr <- effects[chromosome == chr, .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_d)]
  if (!nrow(se_chr)) next

  snplist_chr <- file.path(tmp_dir, sprintf("prs_chr%d.snps", chr))
  fwrite(se_chr[, .(SNP)], snplist_chr, col.names = FALSE)

  bfile      <- sprintf("%s%d_filtered", filtered_prefix, chr)
  out_prefix <- file.path(tmp_dir, sprintf("prs_extract_chr%d", chr))
  cmd <- sprintf(
    "%s --bfile %s --keep %s --extract %s --make-bed --out %s --threads %d",
    shQuote(plink2_path),
    shQuote(bfile),
    shQuote(test_ids_file),
    shQuote(snplist_chr),
    shQuote(out_prefix),
    threads_plink
  )
  ret <- system(cmd)
  if (ret != 0L || !file.exists(paste0(out_prefix, ".bed"))) next

  gd      <- genio::read_plink(out_prefix, verbose = FALSE)
  fam_ids <- gd$fam$id
  snp_ids <- gd$bim$id
  X       <- t(gd$X)

  # Align rows to test IDs
  idx <- match(test_ids_vec, fam_ids)
  X   <- X[idx, , drop = FALSE]

  # Align columns to effects order
  se <- se_chr[match(snp_ids, se_chr$SNP)]
  keep <- !is.na(se$SNP)
  X  <- X[, keep, drop = FALSE]
  se <- se[keep]
  snp_ids <- snp_ids[keep]

  if (ncol(X) == 0) next

  a_eff <- se$a_hat
  d_eff <- se$d_hat

  X <- impute_mode_matrix(X)

  contrib <- numeric(nrow(X))
  for (j in seq_len(ncol(X))) {
    g  <- X[, j]
    a  <- a_eff[j]
    d  <- d_eff[j]
    pd <- se$p_d[j]

    c_add <- numeric(length(g))
    c_add[g == 0L] <- -a
    c_add[g == 1L] <-  0
    c_add[g == 2L] <-  a

    if (!is.na(pd) && pd < dom_p_thresh) {
      c_dom <- numeric(length(g)); c_dom[g == 1L] <- d
      contrib <- contrib + c_add + c_dom
    } else {
      contrib <- contrib + c_add
    }
  }
  prs_vec <- prs_vec + contrib
}

prs_a_d_path <- file.path(pheno_dir, "score", "testing_prs_a_d_method2.txt")
dir.create(dirname(prs_a_d_path), showWarnings = FALSE, recursive = TRUE)
fwrite(data.table(IID = test_ids_vec, PRS = prs_vec), prs_a_d_path, sep = "\t")

ph_test <- fread(test_pheno_file)[, .(IID, pheno)]
prs_dt  <- fread(prs_a_d_path)
res     <- merge(prs_dt, ph_test, by = "IID")
rval    <- suppressWarnings(cor(res$pheno, res$PRS, use = "complete.obs"))

corr_file <- file.path(pheno_dir, "metrics", "corr_ad.txt")
dir.create(dirname(corr_file), showWarnings = FALSE, recursive = TRUE)
write.table(data.frame(corr_ad = rval, N = nrow(res)),
            corr_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
cat("Saved PRS_a_d -> ", prs_a_d_path, "\n", sep = "")
