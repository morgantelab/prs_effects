#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
})

## ---------- Args (must be 9) ----------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)

ids_file          <- args[1]  # FID/IID list
filter_dir        <- args[2]  # directory containing ukb_chr{1..22}_filtered.{bed,bim,fam}
causal_ids_file   <- args[3]  # one-column SNP IDs
true_effects_file <- args[4]  # table with rsid, a, d  (a will be used; d will be overwritten by our rule)
out_pheno         <- args[5]  # output phenotype file
out_geno          <- args[6]  # output genotype-score file
h2                <- as.numeric(args[7])  # heritability
seed              <- as.integer(args[8])  # replicate-specific seed
plink2_path       <- args[9]              # full path to plink2

if (!file.exists(plink2_path)) stop("plink2 not found: ", plink2_path)
if (!dir.exists(filter_dir))   stop("filter_dir missing: ", filter_dir)
if (!file.exists(causal_ids_file)) stop("causal_ids_file missing: ", causal_ids_file)
if (!file.exists(true_effects_file)) stop("true_effects_file missing: ", true_effects_file)

## ---------- Load inputs ----------
ID <- tryCatch(
  fread(ids_file),
  error = function(e) fread(ids_file, header = FALSE, col.names = c("FID","IID"))
)
if (!all(c("FID","IID") %in% names(ID))) {
  setnames(ID, old = names(ID)[1:2], new = c("FID","IID"))
}
ID[, `:=`(FID = as.character(FID), IID = as.character(IID))]

causal_ids <- fread(causal_ids_file, header = FALSE, col.names = "snp")
eff        <- fread(true_effects_file)
stopifnot(all(c("rsid","a","d") %in% names(eff)))
setDT(eff); setkey(eff, rsid)

bed_paths <- sort(list.files(filter_dir, pattern = "^ukb_chr[0-9]+_filtered\\.bed$", full.names = TRUE))
if (length(bed_paths) == 0) stop("No filtered BEDs found in: ", filter_dir)
prefixes <- sub("\\.bed$", "", bed_paths)

## ---------- PASS 1: additive-only scores A_i ----------
# Purpose: compute A_i using a only; later use A_i to stratify individuals into quintiles.
A_dt <- data.frame(FID = ID$FID, IID = ID$IID, A = 0)
set.seed(seed)

for (pref in prefixes) {
  bimf <- paste0(pref, ".bim")
  bedf <- paste0(pref, ".bed")
  famf <- paste0(pref, ".fam")
  if (!file.exists(bimf) || !file.exists(bedf) || !file.exists(famf)) next

  # Intersect current chromosome with causal SNP list
  bim  <- fread(bimf, col.names = c("chr","snp","cm","pos","a1","a2"))
  hits <- intersect(bim$snp, causal_ids$snp)
  if (!length(hits)) next

  snplist_chr <- tempfile(fileext = ".snplist")
  fwrite(data.table(SNP = hits), snplist_chr, col.names = FALSE)

  out_pref <- tempfile(pattern = "extract_")
  cmd <- sprintf(
    "%s --bfile %s --extract %s --make-bed --out %s",
    shQuote(plink2_path), shQuote(pref), shQuote(snplist_chr), shQuote(out_pref)
  )
  ret <- system(cmd)
  if (ret != 0 || !file.exists(paste0(out_pref, ".bed"))) {
    unlink(paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
    next
  }

  fam <- fread(paste0(out_pref, ".fam"),
               col.names = c("FID","IID","pat","mat","sex","pheno"))
  fam[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  idx <- match(fam$IID, A_dt$IID)
  if (any(is.na(idx))) stop("FAM IIDs not all found in ID file for ", out_pref)
  A_dt <- A_dt[idx, , drop = FALSE]
  A_dt$FID <- fam$FID
  A_dt$IID <- fam$IID

  G <- BEDMatrix(paste0(out_pref, ".bed"))
  if (is.null(ncol(G)) || ncol(G) == 0) {
    unlink(paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
    next
  }

  bim2 <- fread(paste0(out_pref, ".bim"),
                col.names = c("chr","snp","pos_cM","pos","alt","ref"))
  eff_block <- eff[J(bim2$snp)]
  eff_block[is.na(a), a := 0]

  for (j in seq_len(ncol(G))) {
    a_j <- eff_block$a[j]
    if (a_j == 0) next
    g <- G[, j]
    contrib <- numeric(length(g))
    contrib[g == 0L] <- -a_j
    contrib[g == 1L] <-  0
    contrib[g == 2L] <-  a_j
    A_dt$A <- A_dt$A + contrib
  }

  unlink(c(snplist_chr, paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex"))), force = TRUE)
}

# Save additive-only genetic value for diagnostics
out_add_only <- sub("\\.([^.]+)$", "_additive_only.txt", out_geno)
dir.create(dirname(out_add_only), recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(FID=A_dt$FID, IID=A_dt$IID, A_additive=A_dt$A),
            out_add_only, row.names = FALSE, quote = FALSE, sep = "\t")

# ---------- Quintiles from additive-only A ----------
# Q1 = lowest additive, Q5 = highest (right tail).
q5 <- cut(rank(A_dt$A, ties.method = "average") / nrow(A_dt),
          breaks = c(0, .2, .4, .6, .8, 1), labels = FALSE, include.lowest = TRUE)
mask_Q1 <- (q5 == 1L)
mask_Q2 <- (q5 == 2L)
mask_Q3 <- (q5 == 3L)
mask_Q4 <- (q5 == 4L)
mask_Q5 <- (q5 == 5L)

## ---------- PASS B: make one fixed d per locus by quintile-specific het presence ----------
# Purpose: decide a single d for each locus based on *presence of heterozygotes* in quintiles.
# Priority: Q5 -> d=+a; else Q1 -> d=-a; else Q2 -> d=-0.8a; else Q4 -> d=+0.8a; else Q3/none -> d=0
RHO <- 1.0  # overall dominance multiplier (kept for convenience if you want to scale)

d_list <- list()

for (pref in prefixes) {
  bimf <- paste0(pref, ".bim"); bedf <- paste0(pref, ".bed"); famf <- paste0(pref, ".fam")
  if (!file.exists(bimf) || !file.exists(bedf) || !file.exists(famf)) next

  # Intersect with causal SNPs for this chromosome
  bim  <- fread(bimf, col.names = c("chr","snp","cm","pos","a1","a2"))
  hits <- intersect(bim$snp, causal_ids$snp)
  if (!length(hits)) next

  # Extract only those SNPs
  snplist_chr <- tempfile(fileext = ".snplist")
  fwrite(data.table(SNP = hits), snplist_chr, col.names = FALSE)
  out_pref <- tempfile(pattern = "extract_")
  cmd <- sprintf('%s --bfile %s --extract %s --make-bed --out %s',
                 shQuote(plink2_path), shQuote(pref), shQuote(snplist_chr), shQuote(out_pref))
  ret <- system(cmd)
  if (ret != 0 || !file.exists(paste0(out_pref, ".bed"))) {
    unlink(c(snplist_chr, paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex"))), force = TRUE)
    next
  }

  # Align masks to this block’s FAM order
  fam <- fread(paste0(out_pref, ".fam"), col.names = c("FID","IID","pat","mat","sex","pheno"))
  fam[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  idx <- match(fam$IID, A_dt$IID)
  if (any(is.na(idx))) stop("FAM IIDs not all found in A_dt for ", out_pref)

  Q1_chr <- mask_Q1[idx]
  Q2_chr <- mask_Q2[idx]
  Q3_chr <- mask_Q3[idx]
  Q4_chr <- mask_Q4[idx]
  Q5_chr <- mask_Q5[idx]

  # Genotypes and effects for this block
  G <- BEDMatrix(paste0(out_pref, ".bed"))
  if (is.null(ncol(G)) || ncol(G) == 0) {
    unlink(c(snplist_chr, paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex"))), force = TRUE)
    next
  }

  bim2 <- fread(paste0(out_pref, ".bim"),
                col.names = c("chr","snp","pos_cM","pos","alt","ref"))
  eff_block <- eff[J(bim2$snp)]
  eff_block[is.na(a), a := 0]

  d_vec <- numeric(ncol(G))

  for (j in seq_len(ncol(G))) {
    a_j <- eff_block$a[j]
    if (a_j == 0) { d_vec[j] <- 0; next }

    g <- G[, j]
    has_Q5 <- any(g[Q5_chr] == 1L, na.rm = TRUE)
    has_Q1 <- any(g[Q1_chr] == 1L, na.rm = TRUE)
    has_Q2 <- any(g[Q2_chr] == 1L, na.rm = TRUE)
    has_Q4 <- any(g[Q4_chr] == 1L, na.rm = TRUE)
    has_Q3 <- any(g[Q3_chr] == 1L, na.rm = TRUE)

    # Priority: Q5 -> Q1 -> Q2 -> Q4 -> Q3/none
    if (has_Q5) {
      d_vec[j] <-  RHO * a_j           # +a
    } else if (has_Q1) {
      d_vec[j] <- -RHO * a_j           # -a
    } else if (has_Q2) {
      d_vec[j] <- -0.8 * RHO * a_j     # -0.8a
    } else if (has_Q4) {
      d_vec[j] <-  0.8 * RHO * a_j     # +0.8a
    } else if (has_Q3) {
      d_vec[j] <-  0                   # 0
    } else {
      d_vec[j] <-  0                   # no hets seen in Q1..Q5
    }
  }

  d_list[[length(d_list) + 1L]] <- data.table(rsid = bim2$snp, a = eff_block$a, d = d_vec)

  unlink(c(snplist_chr, paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex"))), force = TRUE)
}

stopifnot(length(d_list) > 0)
d_all <- rbindlist(d_list, use.names = TRUE, fill = TRUE)
d_all <- d_all[!duplicated(rsid)]
setkey(d_all, rsid)

# Write one fixed d per locus
eff_out <- merge(eff[, .(rsid, a)], d_all[, .(rsid, d)], by = "rsid", all.x = TRUE)
eff_out[is.na(d), d := 0]
out_eff <- sub("\\.([^.]+)$", "_with_d.tsv", true_effects_file)
fwrite(eff_out, out_eff, sep = "\t")
cat("Wrote fixed-d effects to:\n  ", out_eff, "\n", sep = "")

gc()

## ---------- PASS C: simulate genetic value (and phenotype) using fixed d ----------
true_eff <- fread(out_eff)
stopifnot(all(c("rsid","a","d") %in% names(true_eff)))
setDT(true_eff); setkey(true_eff, rsid)

score_dt <- data.frame(FID = ID$FID, IID = ID$IID, score = 0)

for (pref in prefixes) {
  bimf <- paste0(pref, ".bim")
  bedf <- paste0(pref, ".bed")
  famf <- paste0(pref, ".fam")
  if (!file.exists(bimf) || !file.exists(bedf) || !file.exists(famf)) next

  bim  <- fread(bimf, col.names = c("chr","snp","cm","pos","a1","a2"))
  hits <- intersect(bim$snp, causal_ids$snp)
  if (!length(hits)) next

  # write snplist & extract just those SNPs
  snplist_chr <- tempfile(fileext = ".snplist")
  fwrite(data.table(SNP = hits), snplist_chr, col.names = FALSE)

  out_pref <- tempfile(pattern = "extract_")
  cmd <- sprintf(
    "%s --bfile %s --extract %s --make-bed --out %s",
    shQuote(plink2_path),
    shQuote(pref),
    shQuote(snplist_chr),
    shQuote(out_pref)
  )
  ret <- system(cmd)
  if (ret != 0 || !file.exists(paste0(out_pref, ".bed"))) {
    # clean and skip if extract failed/empty
    unlink(paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
    next
  }

  fam <- fread(paste0(out_pref, ".fam"),
               col.names = c("FID","IID","pat","mat","sex","pheno"))
  fam[, `:=`(FID = as.character(FID), IID = as.character(IID))]

  # Reorder score_dt to FAM order (robust alignment)
  idx <- match(fam$IID, score_dt$IID)
  if (any(is.na(idx))) stop("FAM IIDs not all found in ID file for ", out_pref)
  score_dt <- score_dt[idx, , drop = FALSE]
  # overwrite FID/IID to exactly match FAM (safety)
  score_dt$FID <- fam$FID
  score_dt$IID <- fam$IID

  G <- BEDMatrix(paste0(out_pref, ".bed"))   # N_ind x N_snp
  if (is.null(ncol(G)) || ncol(G) == 0) {
    unlink(paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex")), force = TRUE)
    next
  }

  bim2 <- fread(paste0(out_pref, ".bim"),
                col.names = c("chr","snp","pos_cM","pos","alt","ref"))

  # Vectorized look-up of effects for all SNPs in this extract
  eff_block <- true_eff[J(bim2$snp)]
  eff_block[is.na(a), a := 0]
  eff_block[is.na(d), d := 0]

  # Accumulate contribution for each SNP column
  # mapping: g=0 -> -a; g=1 -> d; g=2 -> +a
  for (j in seq_len(ncol(G))) {
    a <- eff_block$a[j]
    d <- eff_block$d[j]
    if ((a == 0) && (d == 0)) next
    g <- G[, j]
    contrib <- numeric(length(g))
    contrib[g == 0L] <- -a
    contrib[g == 1L] <-  d
    contrib[g == 2L] <-  a
    score_dt$score <- score_dt$score + contrib
  }

  # clean temp files for this chromosome
  unlink(c(snplist_chr,
           paste0(out_pref, c(".bed",".bim",".fam",".log",".nosex"))), force = TRUE)
}

## ---------- Phenotype ----------
var_g <- var(score_dt$score)
if (!is.finite(var_g) || var_g <= 0) stop("Genetic score variance is non-positive.")
var_e <- ((1 / h2) - 1) * var_g

set.seed(seed)
pheno_vec <- if (is.finite(var_e) && var_e > 0) {
  score_dt$score + rnorm(nrow(score_dt), mean = 0, sd = sqrt(var_e))
} else {
  score_dt$score
}

ph <- data.frame(FID = as.character(score_dt$FID),
                 IID = as.character(score_dt$IID),
                 pheno = pheno_vec)

dir.create(dirname(out_pheno), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_geno),  recursive = TRUE, showWarnings = FALSE)

write.table(ph, out_pheno, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(score_dt[, c("FID","IID","score")], out_geno, row.names = FALSE, quote = FALSE, sep = "\t")
