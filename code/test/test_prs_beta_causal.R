##-------------------------------------------
## ---- PRS using ONLY causal loci (no clumping) ----
# Inputs already created above:
#   - output_causal_snps_IDs : list of causal rsIDs (one per line)
#   - true_effects_file      : data.table with columns rsid, a, d
#   - pheno_dir              : base working directory with ukb_chr{1..22}_filtered.{bed/bim/fam}
#   - output_testing_ID      : testing IDs (IID/FID two columns, no header)
#   - output_testing_pheno   : testing phenotype table with FID,IID,pheno

# Load causal IDs and true effects
causal_ids_dt <- fread(output_causal_snps_IDs, header = FALSE, col.names = "ID")
#-------------------------------
# (B) GWAS β restricted to causal rsIDs (no clumping)
#   Upper bound that still uses estimated effects; weights = training GWAS betas on causal rsIDs only.
#-------------------------------
# --- Causal-only PRS using GWAS betas: train-on-train, predict-on-test ---

# 1) 权重：仅限因果 rsIDs
gwas_all <- fread(file.path(pheno_dir, "gwas_merged.linear"))
if ("#CHROM" %in% names(gwas_all)) setnames(gwas_all, "#CHROM", "CHR")
if ("TEST" %in% names(gwas_all))   gwas_all <- gwas_all[TEST == "ADD"]

causal_ids_dt <- fread(output_causal_snps_IDs, header = FALSE, col.names = "ID")
gwas_causal   <- merge(causal_ids_dt, gwas_all[, .(ID, A1, BETA)], by = "ID")
wgt_gwas_path <- file.path(pheno_dir, "beta_for_score.causal_gwas.txt")
fwrite(gwas_causal, wgt_gwas_path, sep = "\t")


# 2) 测试集逐染色体打分并合并

for (chr in 1:22) {
  bfile <- file.path(pheno_dir, sprintf("ukb_chr%d_filtered", chr))
  if (!file.exists(paste0(bfile, ".bed"))) next
  cmd <- sprintf(
    "%s --bfile %s --keep %s --score %s 1 2 3 header-read cols=+scoresums --out %s",
    shQuote(plink2_path),
    shQuote(bfile),
    shQuote(output_testing_ID),
    shQuote(wgt_gwas_path),
    shQuote(file.path(pheno_dir, sprintf("score_causalGWAS_test_chr%d", chr)))
  )
  system(cmd)
}

ss_test <- file.path(pheno_dir, sprintf("score_causalGWAS_test_chr%d.sscore", 1:22))
ss_test <- ss_test[file.exists(ss_test)]
stopifnot(length(ss_test) > 0)
dt0_te  <- fread(ss_test[1])
col_sum <- "BETA_SUM"
prs_te  <- dt0_te[[col_sum]]
for (f in ss_test[-1]) prs_te <- prs_te + fread(f)[[col_sum]]

test_prs_dt <- data.table(IID = dt0_te$IID, PRS = prs_te)
fwrite(test_prs_dt, file.path(pheno_dir, "testing_prs_causal_gwas.txt"), sep = "\t")

# 5) 用训练得到的模型，在测试集预测并评估
test_pheno <- fread(output_testing_pheno)[, .(IID, pheno)]
test_prs_dt <- fread(file.path(pheno_dir, "testing_prs_causal_gwas.txt"))
test_m     <- merge(test_prs_dt, test_pheno, by = "IID")

r_te   <- suppressWarnings(cor(test_m$pheno, test_m$PRS, use = "complete.obs"))
r_te

#> r_te
#0.6526649
# > r_te
# [1] 0.6832614
# > r_te
# [1] 0.6986369

#100comdom
# > r_te
# [1] 0.627736
# 50com50add
# > r_te
# [1] 0.6825314
# 100 add
# > r_te
# [1] 0.706936