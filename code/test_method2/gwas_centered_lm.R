#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
  library(genio)
})

# -------------------- Args --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: gwas_centered_lm.R ",
       "FILTERED_PREFIX CHR TRAIN_IDS TRAIN_PHENO PHENO_NAME OUT_TSV CHUNK_SIZE")
}

filtered_prefix <- args[1]   # e.g. /.../OUT_DIR/filtered/ukb_chr
chr             <- as.integer(args[2])  # e.g. 1..22
train_ids_file  <- args[3]
train_pheno_file<- args[4]
pheno_name      <- args[5]
out_tsv         <- args[6]
chunk_size      <- as.integer(args[7])

# -------------------- Helpers --------------------
fail <- function(...) { message(...); quit(status = 1) }

# Input files
bed_pref <- sprintf("%s%d_filtered", filtered_prefix, chr)
bed_fp <- paste0(bed_pref, ".bed")
bim_fp <- paste0(bed_pref, ".bim")
fam_fp <- paste0(bed_pref, ".fam")
for (f in c(bed_fp, bim_fp, fam_fp)) if (!file.exists(f)) fail("Missing file: ", f)

# Read BIM/FAM
bim <- fread(bim_fp, header = FALSE)
setnames(bim, c("CHR","ID","CM","POS","A1","A2"))
fam <- fread(fam_fp, header = FALSE)
setnames(fam, c("FID","IID","PID","MID","SEX","PHENO"))

# Training IDs and phenotype
train_ids <- fread(train_ids_file)
if (ncol(train_ids) < 2) fail("TRAIN_IDS must have at least two columns: FID IID")
setnames(train_ids, c("FID","IID", names(train_ids)[-(1:2)]))

train_ph <- fread(train_pheno_file)
if (!(pheno_name %in% names(train_ph))) fail("Column not found in TRAIN_PHENO: ", pheno_name)

# Align FAM and phenotype
fam[, ROW := .I]
setkey(fam, FID, IID)
setkey(train_ids, FID, IID)
tr <- fam[train_ids, nomatch = 0L]
if (nrow(tr) == 0) fail("No overlap between FAM and TRAIN_IDS.")

setkey(train_ph, FID, IID)
tr_ph <- train_ph[tr[, .(FID, IID)], nomatch = 0L]
if (nrow(tr_ph) != nrow(tr)) fail("TRAIN_PHENO does not fully cover TRAIN_IDS.")
y <- as.numeric(tr_ph[[pheno_name]])
if (any(!is.finite(y))) fail("Non-finite phenotype values found.")

# Load genotype (BEDMatrix)
G <- BEDMatrix(bed_fp)
G <- G[tr$ROW, , drop = FALSE]  # subset to training individuals

n <- nrow(G)
m <- ncol(G)
if (m == 0) fail("No SNPs for this chromosome.")

message(sprintf("Running lm-based GWAS on chr %d: n=%d, m=%d", chr, n, m))

# Prepare output table
res <- data.table(
  `#CHROM` = as.integer(bim$CHR),
  POS      = as.integer(bim$POS),
  ID       = bim$ID,
  REF      = bim$A2,  # treat A2 as reference
  ALT      = bim$A1,  # treat A1 as effect allele
  A1       = bim$A1,
  TEST     = rep("ADD", m),
  OBS_CT   = integer(m),
  BETA     = rep(NA_real_, m),
  SE       = rep(NA_real_, m),
  T_STAT   = rep(NA_real_, m),
  P        = rep(NA_real_, m)
)

# Function for one SNP
fit_one <- function(x, y) {
  ok <- is.finite(y) & !is.na(x)
  n_ok <- sum(ok)
  if (n_ok < 3L) return(list(n=n_ok, beta=NA_real_, se=NA_real_, t=NA_real_, p=NA_real_))
  x_c <- x[ok] - mean(x[ok])
  y_s <- y[ok]
  fit <- try(lm(y_s ~ x_c), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(n=n_ok, beta=NA_real_, se=NA_real_, t=NA_real_, p=NA_real_))
  co <- summary(fit)$coefficients
  if (nrow(co) < 2) return(list(n=n_ok, beta=NA_real_, se=NA_real_, t=NA_real_, p=NA_real_))
  list(
    n    = n_ok,
    beta = unname(co[2, "Estimate"]),
    se   = unname(co[2, "Std. Error"]),
    t    = unname(co[2, "t value"]),
    p    = unname(co[2, "Pr(>|t|)"])
  )
}

# Main loop (chunked)
start <- 1L
while (start <= m) {
  end <- min(start + chunk_size - 1L, m)
  cols <- start:end
  X <- G[, cols, drop = FALSE]
  for (j in seq_along(cols)) {
    idx <- cols[j]
    r <- fit_one(as.numeric(X[, j]), y)
    res$OBS_CT[idx] <- r$n
    res$BETA[idx]   <- r$beta
    res$SE[idx]     <- r$se
    res$T_STAT[idx] <- r$t
    res$P[idx]      <- r$p
  }
  message(sprintf(".. processed %d / %d SNPs", end, m))
  start <- end + 1L
}

# Write output
dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
fwrite(res, out_tsv, sep = "\t", quote = FALSE, na = "NA")
message("Completed writing: ", out_tsv)
