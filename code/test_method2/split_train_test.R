#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly=TRUE)
ids_file     <- args[1]
pheno_file   <- args[2]
out_train_ids<- args[3]
out_test_ids <- args[4]
out_train_ph <- args[5]
out_test_ph  <- args[6]
SEED         <- as.integer(args[7])

ID <- fread(ids_file, header=FALSE, col.names=c("FID","IID"))
ph <- fread(pheno_file)

set.seed(SEED)
perm <- sample(nrow(ID))
train <- ID[perm[1:180000]]
test  <- ID[perm[180001:200000]]

fwrite(train, out_train_ids, sep="\t", col.names=FALSE)
fwrite(test,  out_test_ids,  sep="\t", col.names=FALSE)

setDT(ph)
fwrite(ph[IID %in% train$IID], out_train_ph, sep="\t")
fwrite(ph[IID %in% test$IID],  out_test_ph,  sep="\t")
