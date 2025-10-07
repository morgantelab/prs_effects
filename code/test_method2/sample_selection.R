#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})
args <- commandArgs(trailingOnly=TRUE)
fam_file    <- args[1]
pheno_file  <- args[2]
out_ids     <- args[3]
N_SAMPLE    <- as.integer(args[4])
SAMPLE_SEED <- as.integer(args[5])

fam <- fread(fam_file, header=FALSE,
             col.names=c("FID","IID","FatherID","MotherID","Sex","Phenotype"))
keep_ids <- fam[Sex %in% c(1,2), .(IID)]

ph  <- fread(pheno_file)
phf <- ph %>% select(
  eid,
  age = `21022-0.0`,
  sex = `31-0.0`,
  pc1 = `22009-0.1`, pc2 = `22009-0.2`, pc3 = `22009-0.3`, pc4 = `22009-0.4`, pc5 = `22009-0.5`,
  pc6 = `22009-0.6`, pc7 = `22009-0.7`, pc8 = `22009-0.8`, pc9 = `22009-0.9`, pc10 = `22009-0.10`,
  used_in_pca = `22020-0.0`, outlier = `22027-0.0`, ethnic = `22006-0.0`
) %>% filter(ethnic==1, used_in_pca==1, eid %in% keep_ids$IID)

set.seed(SAMPLE_SEED)
idx <- sample(nrow(phf), N_SAMPLE)
ID  <- data.frame(FID=phf$eid[idx]); ID$IID <- ID$FID
fwrite(ID, out_ids, sep="\t", col.names=FALSE)
cat("Unique IDs in final file:", nrow(ID), "\n")
