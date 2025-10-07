#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

## ---- Paths / config ----
rep <- 2
results_root <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4"
pheno_file   <- file.path(results_root, "pheno", sprintf("simulated_pheno_rep%d.txt", rep))
out_dir      <- file.path(results_root, "test")

prs_score_path       <- file.path(out_dir, sprintf("test_prs_score_rep%d.txt", rep))
prs_scan_path        <- file.path(out_dir, sprintf("test_prs_scan_rep%d.txt", rep))
prs_gwas_path        <- file.path(out_dir, sprintf("test_prs_gwas_rep%d.txt", rep))
prs_scan_impute_path <- file.path(out_dir, sprintf("test_prs_scan_impute_rep%d.txt", rep))

merged_out     <- file.path(out_dir, sprintf("merged_pheno_prs_rep%d.txt", rep))

## ---- helpers ----
norm_id <- function(x) as.character(trimws(x))   # 统一去空格并转字符
first_prs_col <- function(dt) setdiff(names(dt), c("FID","IID"))[1]

## ---- Load phenotype ----
pheno <- fread(pheno_file)                          # 期望包含 FID, IID, pheno
pheno[, `:=`(FID = norm_id(FID), IID = norm_id(IID))]

## ---- Load PRS_score ----
prs_score <- fread(prs_score_path)
prs_score[, IID := norm_id(IID)]
setnames(prs_score, first_prs_col(prs_score), "PRS_score")
prs_score <- prs_score[, .(IID, PRS_score)]

## ---- Load PRS_scan ----
prs_scan <- fread(prs_scan_path)
prs_scan[, IID := norm_id(IID)]
setnames(prs_scan, first_prs_col(prs_scan), "PRS_scan")
prs_scan <- prs_scan[, .(IID, PRS_scan)]

## ---- Load PRS_scan_impute ----
prs_scan_impute <- fread(prs_scan_impute_path)
prs_scan_impute[, IID := norm_id(IID)]
setnames(prs_scan_impute, first_prs_col(prs_scan_impute), "PRS_scan_impute")
prs_scan_impute <- prs_scan_impute[, .(IID, PRS_scan_impute)]

## ---- Load PRS_gwas ----
prs_gwas <- fread(prs_gwas_path)
prs_gwas[, IID := norm_id(IID)]
setnames(prs_gwas, first_prs_col(prs_gwas), "PRS_gwas")
prs_gwas <- prs_gwas[, .(IID, PRS_gwas)]

## ---- Quick overlap diagnostics ----
n_pheno <- nrow(pheno)
ov <- function(x) length(intersect(pheno$IID, x$IID))
cat(sprintf("Samples in pheno: %d\n", n_pheno))
cat(sprintf("Overlap with PRS_scan:         %d\n", ov(prs_scan)))
cat(sprintf("Overlap with PRS_scan_impute:  %d\n", ov(prs_scan_impute)))
cat(sprintf("Overlap with PRS_gwas:         %d\n", ov(prs_gwas)))
cat(sprintf("Overlap with PRS_score:        %d\n", ov(prs_score)))


## ---- Merge by IID (left join to pheno) ----
merged <- Reduce(function(x, y) merge(x, y, by = "IID"),
                 list(pheno[, .(FID, IID, pheno)],
                      prs_scan, prs_scan_impute, prs_gwas, prs_score))

## ---- Order columns nicely ----
setcolorder(merged, c("FID","IID","pheno",
                      "PRS_scan","PRS_scan_impute","PRS_scan_final",
                      "PRS_gwas","PRS_score"))

## ---- Save ----
fwrite(merged, merged_out, sep = "\t")
cat("✅ Saved:", merged_out, "\n")

## ---- Quick sanity (optional) ----
avail <- names(merged)
if ("PRS_scan" %in% avail)         cat("cor(pheno, PRS_scan)         =", cor(merged$pheno, merged$PRS_scan,         use="complete.obs"), "\n")
if ("PRS_scan_impute" %in% avail)  cat("cor(pheno, PRS_scan_impute)  =", cor(merged$pheno, merged$PRS_scan_impute,  use="complete.obs"), "\n")
if ("PRS_gwas" %in% avail)         cat("cor(pheno, PRS_gwas)         =", cor(merged$pheno, merged$PRS_gwas,         use="complete.obs"), "\n")
if ("PRS_score" %in% avail)        cat("cor(pheno, PRS_score)        =", cor(merged$pheno, merged$PRS_score,        use="complete.obs"), "\n")
