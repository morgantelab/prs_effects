#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(genio) })
args <- commandArgs(trailingOnly=TRUE)
effects_file <- args[1]
test_ids_file<- args[2]
test_ph_file <- args[3]
out_prs      <- args[4]
out_corr     <- args[5]
FILTER_DIR   <- args[6]
DOM_PTHRESH  <- as.numeric(args[7])

effects <- fread(effects_file)
effects[, chromosome := as.integer(chromosome)]
ids <- fread(test_ids_file, header=FALSE)[[2]]
ph  <- fread(test_ph_file)[, .(IID=as.character(IID), pheno)]

impute_mode_matrix <- function(X) {
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[,j])) {
      g <- X[,j]
      tab <- tabulate(g[!is.na(g)] + 1L, nbins=3L)
      g[is.na(g)] <- which.max(tab) - 1L
      X[,j] <- g
    }
  }
  X
}

prs <- numeric(length(ids)); names(prs) <- ids

for (chr in sort(unique(effects$chromosome))) {
  se <- effects[chromosome==chr, .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_d)]
  if (!nrow(se)) next

  snpfile <- tempfile(fileext=".snps")
  fwrite(se[, .(SNP)], snpfile, col.names=FALSE)

  pref <- file.path(FILTER_DIR, sprintf("ukb_chr%d_filtered", chr))
  outp <- tempfile(pattern="prs_extract_")
  cmd <- sprintf("plink2 --bfile %s --keep %s --extract %s --make-bed --out %s",
                 pref, test_ids_file, snpfile, outp)
  system(cmd)

  if (!file.exists(paste0(outp, ".bed"))) next

  gd <- genio::read_plink(outp, verbose=FALSE)
  fam_ids <- gd$fam$id
  snp_ids <- gd$bim$id
  X <- t(gd$X)  # N_ind x N_snp

  idx <- match(ids, fam_ids)
  X <- X[idx, , drop=FALSE]

  se2 <- se[match(snp_ids, se$SNP)]
  keep <- !is.na(se2$SNP)
  X <- X[, keep, drop=FALSE]
  se2 <- se2[keep]

  a  <- se2$a_hat
  d  <- se2$d_hat
  pd <- se2$p_d

  X <- impute_mode_matrix(X)

  contrib <- numeric(nrow(X))
  for (j in seq_len(ncol(X))) {
    g <- X[, j]
    aa <- a[j]
    dd <- d[j]
    p  <- pd[j]
    add <- ifelse(g==0, -aa, ifelse(g==1, 0, ifelse(g==2, aa, 0)))
    if (!is.na(p) && p < DOM_PTHRESH) add <- add + ifelse(g==1, dd, 0)
    contrib <- contrib + add
  }
  prs <- prs + contrib
}

prs_dt <- data.table(IID=ids, PRS=prs)
fwrite(prs_dt, out_prs, sep="\t")
res <- merge(prs_dt, ph, by="IID")
rval <- suppressWarnings(cor(res$pheno, res$PRS, use="complete.obs"))
write.table(data.frame(corr_ad=rval, N=nrow(res)),
            out_corr, row.names=FALSE, quote=FALSE, sep="\t")
cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
