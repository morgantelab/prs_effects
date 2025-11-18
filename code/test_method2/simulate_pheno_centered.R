#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
})

## ---------- Args (must be 9) ----------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 9)

ids_file          <- args[1]  # FID/IID list (likely with header)
filter_dir        <- args[2]  # directory containing ukb_chr{1..22}_filtered.{bed,bim,fam}
causal_ids_file   <- args[3]  # one-column SNP IDs
true_effects_file <- args[4]  # table with rsid, a, d
out_pheno         <- args[5]  # output phenotype file
out_geno          <- args[6]  # output genotype-score file
h2                <- as.numeric(args[7])  # heritability
seed              <- as.integer(args[8])  # replicate-specific seed
plink2_path       <- args[9]              # full path to plink2

if (!file.exists(plink2_path)) stop("plink2 not found: ", plink2_path)
if (!dir.exists(filter_dir))   stop("filter_dir missing: ", filter_dir)
if (!file.exists(causal_ids_file)) stop("causal_ids_file missing: ", causal_ids_file)
if (!file.exists(true_effects_file)) stop("true_effects_file missing: ", true_effects_file)

## ---------- Inputs ----------
# IDs: try headered, else fallback
ID <- tryCatch(
  fread(ids_file),
  error = function(e) fread(ids_file, header = FALSE, col.names = c("FID","IID"))
)
if (!all(c("FID","IID") %in% names(ID))) {
  # fallback if different header case
  setnames(ID, old = names(ID)[1:2], new = c("FID","IID"))
}
ID[, `:=`(FID = as.character(FID), IID = as.character(IID))]

causal_ids <- fread(causal_ids_file, header = FALSE, col.names = "snp")
true_eff   <- fread(true_effects_file)
stopifnot(all(c("rsid","a","d") %in% names(true_eff)))
setDT(true_eff); setkey(true_eff, rsid)

## ---------- Collect per-chr prefixes from directory ----------
bed_paths <- sort(list.files(filter_dir, pattern = "^ukb_chr[0-9]+_filtered\\.bed$", full.names = TRUE))
if (length(bed_paths) == 0) stop("No filtered BEDs found in: ", filter_dir)
prefixes <- sub("\\.bed$", "", bed_paths)  # "/.../ukb_chr{c}_filtered"

## ---------- Init ----------
score_dt <- data.frame(FID = ID$FID, IID = ID$IID, score = 0)
set.seed(seed)

## ---------- Loop chromosomes present in filter_dir ----------
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
  # requires key 'rsid' and vector of SNP IDs
  eff_block <- true_eff[J(bim2$snp)]
  # If any missing, set a/d to 0 so they don't contribute
  eff_block[is.na(a), a := 0]
  eff_block[is.na(d), d := 0]

  # Accumulate contribution for each SNP column
  # mapping: g=0 -> -a; g=1 -> d; g=2 -> +a
  for (j in seq_len(ncol(G))) {
    a <- eff_block$a[j]
    d <- eff_block$d[j]
    if ((a == 0) && (d == 0)) next
    g <- G[, j]
    # fast mapping without branching
    # contrib = (-a)*(g==0) + d*(g==1) + (+a)*(g==2)
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

## ---------- Simulate phenotype ----------
var_g <- var(score_dt$score)
if (!is.finite(var_g) || var_g <= 0) stop("Genetic score variance is non-positive.")
var_e <- ((1 / h2) - 1) * var_g
set.seed(seed)
score_dt$pheno <- score_dt$score + rnorm(nrow(score_dt), mean = 0, sd = sqrt(var_e))

## ---------- Center phenotype (mean 0) ----------
ph_mean <- mean(score_dt$pheno)
score_dt$pheno_centered <- score_dt$pheno - ph_mean

## ---------- Write outputs ----------
ph_raw <- data.frame(FID = as.character(score_dt$FID),
                     IID = as.character(score_dt$IID),
                     pheno = score_dt$pheno)

ph_center <- data.frame(FID = as.character(score_dt$FID),
                        IID = as.character(score_dt$IID),
                        pheno = score_dt$pheno_centered)

# 1) original (unchanged) phenotype at out_pheno_raw
out_pheno_raw <- sub("(\\.txt)?$", "_raw.txt", out_pheno)

write.table(ph_raw, out_pheno_raw, row.names = FALSE, quote = FALSE, sep = "\t")
# 2) centered phenotype at out_pheno
write.table(ph_center,  out_pheno, row.names = FALSE, quote = FALSE, sep = "\t")

# genotype score
write.table(score_dt[, c("FID","IID","score")], out_geno, row.names = FALSE, quote = FALSE, sep = "\t")
