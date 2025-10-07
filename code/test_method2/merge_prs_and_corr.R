#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly=TRUE)
score_dir <- args[1]
out_prs   <- args[2]
test_ph   <- args[3]
out_corr  <- args[4]

files <- list.files(score_dir, pattern="score_clumped_chr.*\\.sscore$", full.names=TRUE)
stopifnot(length(files) > 0)
dt0 <- fread(files[1])
sum_col <- "BETA_SUM"
prs <- dt0[[sum_col]]
if (length(files) > 1) for (f in files[-1]) prs <- prs + fread(f)[[sum_col]]

prs_dt <- data.table(IID=dt0$IID, PRS=prs)
fwrite(prs_dt, out_prs, sep="\t")

ph <- fread(test_ph)[, .(IID, pheno)]
res <- merge(prs_dt, ph, by="IID")
rval <- suppressWarnings(cor(res$pheno, res$PRS, use="complete.obs"))

write.table(data.frame(corr_beta=rval, N=nrow(res)),
            out_corr, row.names=FALSE, quote=FALSE, sep="\t")
cat(sprintf("[TEST clumped β] N=%d  cor=%.6f\n", nrow(res), rval))
