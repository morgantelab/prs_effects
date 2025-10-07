## ---- GWAS on training ----
pheno_name <- "pheno"
threads_plink <- 20 

for (CHR in 1:22) {
  cmd <- sprintf(
    "%s --bfile %s --pheno %s --keep %s --glm no-x-sex allow-no-covars omit-ref --no-input-missing-phenotype --pheno-name %s --out %s --threads %d",
    shQuote(plink2_path),
    shQuote(file.path(pheno_dir, sprintf("ukb_chr%d_filtered", CHR))),
    shQuote(output_training_pheno),
    shQuote(output_training_ID),
    pheno_name,
    shQuote(file.path(pheno_dir, sprintf("gwas_chr%d", CHR))),
    threads_plink
  )
  system(cmd)
}

gwas_files <- file.path(pheno_dir, sprintf("gwas_chr%d.%s.glm.linear", 1:22, pheno_name))
gwas_files <- gwas_files[file.exists(gwas_files)]
gwas <- data.table::rbindlist(lapply(gwas_files, data.table::fread), use.names = TRUE, fill = TRUE)
data.table::fwrite(gwas, file.path(pheno_dir, "gwas_merged.linear"), sep = "\t")
## ---- Clumping (R2=0.3, KB=250, P1=1e-4) ----
clump_p1 <- 1e-4
clump_r2 <- 0.3
clump_kb <- 250

gwas <- data.table::fread(file.path(pheno_dir, "gwas_merged.linear"))
if ("#CHROM" %in% names(gwas)) data.table::setnames(gwas, "#CHROM", "CHR")

for (chr in 1:22) {
  g_chr <- gwas[CHR == chr, .(ID, P)]
  sumfile   <- file.path(pheno_dir, sprintf("gwas_for_clump_chr%d.txt", chr))
  outprefix <- file.path(pheno_dir, sprintf("clump_chr%d", chr))
  data.table::fwrite(g_chr, sumfile, sep = "\t")
  
  cmd <- sprintf(
    "%s --bfile %s --keep %s --clump %s --clump-snp-field ID --clump-field P --clump-p1 %g --clump-r2 %g --clump-kb %d --out %s",
    shQuote(plink2_path),
    shQuote(file.path(pheno_dir, sprintf("ukb_chr%d_filtered", chr))),
    shQuote(output_training_ID),
    shQuote(sumfile),
    clump_p1, clump_r2, clump_kb,
    shQuote(outprefix)
  )
  system(cmd)
}


clumped_files <- file.path(pheno_dir, sprintf("clump_chr%d.clumps", 1:22))
clump_all <- data.table::rbindlist(lapply(clumped_files[file.exists(clumped_files)],
                                          data.table::fread, fill = TRUE),
                                   use.names = TRUE, fill = TRUE)

idx_snps <- unique(clump_all$ID)
snplist_path <- file.path(pheno_dir, "clumped_snps.snplist")
data.table::fwrite(data.table::data.table(SNP = idx_snps),
                   snplist_path, sep = "\t", col.names = FALSE)
cat("Clumped index SNPs → ", snplist_path, " (N=", length(idx_snps), ")\n", sep = "")

## ---- Use clumped SNPs to calculate PRS in Testing Set ----
# 1) 构建权重：只保留 clumped SNPs（ID/A1/BETA）
gwas_all <- data.table::fread(file.path(pheno_dir, "gwas_merged.linear"))
if ("#CHROM" %in% names(gwas_all)) data.table::setnames(gwas_all, "#CHROM", "CHR")
if ("TEST" %in% names(gwas_all))   gwas_all <- gwas_all[TEST == "ADD"]

idx_snps <- data.table::fread(file.path(pheno_dir, "clumped_snps.snplist"),
                              header = FALSE)[[1]]
wgt_clumped <- gwas_all[ID %in% idx_snps, .(ID, A1, BETA)]
wgt_clumped_path <- file.path(pheno_dir, "beta_for_score.clumped.txt")
data.table::fwrite(wgt_clumped, wgt_clumped_path, sep = "\t")

# 2) 逐染色体打分（PLINK2 会自动忽略当前 bfile 中不存在的 SNP）
for (CHR in 1:22) {
  cmd <- sprintf(
    "%s --bfile %s --keep %s --score %s 1 2 3 header-read cols=+scoresums --out %s",
    shQuote(plink2_path),
    shQuote(file.path(pheno_dir, sprintf("ukb_chr%d_filtered", CHR))),
    shQuote(output_testing_ID),
    shQuote(wgt_clumped_path),
    shQuote(file.path(pheno_dir, sprintf("score_clumped_chr%d", CHR)))
  )
  system(cmd)
}

# 3) 合并 22 条 .sscore → 总 PRS（列名通常为 BETA_SUM）
ss_files <- file.path(pheno_dir, sprintf("score_clumped_chr%d.sscore", 1:22))
ss_files <- ss_files[file.exists(ss_files)]

dt0     <- data.table::fread(ss_files[1])
sum_col <- "BETA_SUM"
prs     <- dt0[[sum_col]]
for (f in ss_files[-1]) prs <- prs + data.table::fread(f)[[sum_col]]

prs_dt <- data.table::data.table(IID = dt0$IID, PRS = prs)
data.table::fwrite(prs_dt, file.path(pheno_dir, "testing_prs_beta_clumped.txt"), sep = "\t")

# 4) （可选）与 testing 表型求相关
testing_pheno <- fread(output_testing_pheno)
testing_pheno$IID <- as.numeric(testing_pheno$IID)
prs_dt <- fread(file.path(pheno_dir, "testing_prs_beta_clumped.txt"))
res     <- merge(prs_dt, testing_pheno, by = "IID")
rval    <- suppressWarnings(cor(res$pheno, res$PRS, use = "complete.obs"))
cat(sprintf("[TEST clumped β] N=%d  cor=%.6f\n", nrow(res), rval))
#[TEST clumped β] N=20000  cor=0.545346
#[TEST clumped β] N=20000  cor=0.547299
#[TEST clumped β] N=20000  cor=0.555227
# 100 comdom
# > cat(sprintf("[TEST clumped β] N=%d  cor=%.6f\n", nrow(res), rval))
# [TEST clumped β] N=20000  cor=0.484090
#50com50add
# > cat(sprintf("[TEST clumped β] N=%d  cor=%.6f\n", nrow(res), rval))
# [TEST clumped β] N=20000  cor=0.527883
#100add
# [TEST clumped β] N=20000  cor=0.544260