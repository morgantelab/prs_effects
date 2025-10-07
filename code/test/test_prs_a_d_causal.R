
##-----------------------------------------
## --- PRS_a_d on ALL causal SNPs (NO clumping, NO allele alignment) ---
library(data.table)
library(genio)

# 依赖对象/路径（与您现有脚本一致）：
# - effects_out                : 训练集效应表 (含 SNP, CHR, ref_allele, alt_allele, a_hat, d_hat, p_d)
# - output_causal_snps_IDs     : 所有因果 rsID（单列）
# - geno_pattern               : sprintf("%s%d_filtered", geno_pattern, chr)
# - output_testing_ID          : 测试集 keep 文件（两列 FID IID，无表头）
# - output_testing_pheno       : 测试集表型 (含 IID, pheno)
# - plink2_path, pheno_dir

# 读入因果位点与训练效应
causal_ids_dt <- fread(output_causal_snps_IDs, header = FALSE, col.names = "SNP")

effects_out <- file.path(pheno_dir, "effects_training.txt")
eff_all       <- fread(effects_out)
eff_all[, CHR := as.integer(chromosome)]

# 仅保留因果位点，且有 a_hat / d_hat
eff_causal <- eff_all[
  SNP %in% causal_ids_dt$SNP & !is.na(a_hat) & !is.na(d_hat),
  .(SNP, CHR, a_hat, d_hat, p_d)
]

# 测试集 IID 顺序
test_ids <- fread(output_testing_ID, header = FALSE)[[2]]

# 众数填补函数
impute_mode_matrix <- function(X) {
  nc <- ncol(X)
  for (j in seq_len(nc)) {
    if (anyNA(X[, j])) {
      g <- X[, j]
      tab <- tabulate(g[!is.na(g)] + 1L, nbins = 3L)
      g[is.na(g)] <- which.max(tab) - 1L
      X[, j] <- g
    }
  }
  X
}

dominance_p_threshold <- 10E-5
prs_vec <- numeric(length(test_ids)); names(prs_vec) <- test_ids
tmp_dir <- file.path(pheno_dir, "tmp_prs_ad_noclump_noalign")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

for (chr in 1:22) {
  se_chr <- eff_causal[CHR == chr, .(SNP, a_hat, d_hat, p_d)]
  if (!nrow(se_chr)) next
  
  # 该染色体的因果 SNP 列表
  snplist_chr <- file.path(tmp_dir, sprintf("prs_noclump_chr%d.snps", chr))
  fwrite(se_chr[, .(SNP)], snplist_chr, col.names = FALSE)
  
  # 提取测试集基因型
  bfile      <- sprintf("%s%d_filtered", geno_pattern, chr)
  out_prefix <- file.path(tmp_dir, sprintf("prs_noclump_extract_chr%d", chr))
  cmd <- sprintf("%s --bfile %s --keep %s --extract %s --make-bed --out %s",
                 shQuote(plink2_path),
                 shQuote(bfile),
                 shQuote(output_testing_ID),
                 shQuote(snplist_chr),
                 shQuote(out_prefix))
  system(cmd)
  if (!file.exists(paste0(out_prefix, ".bed"))) next
  
  gd      <- genio::read_plink(out_prefix, verbose = FALSE)
  fam_ids <- gd$fam$id
  snp_ids <- gd$bim$id
  X       <- t(gd$X)  # N_ind x N_snp（直接按 genio 的 A1 计数，不做等位基因检查/翻转）
  
  # 行顺序对齐
  idx <- match(test_ids, fam_ids)
  X   <- X[idx, , drop = FALSE]
  
  # 列顺序对齐（按 SNP ID 匹配；不检查 A1/A2）
  se <- se_chr[match(snp_ids, se_chr$SNP)]
  keep <- !is.na(se$SNP)
  X       <- X[, keep, drop = FALSE]
  se      <- se[keep]
  snp_ids <- snp_ids[keep]
  
  a_eff <- se$a_hat; a_eff[is.na(a_eff)] <- 0
  d_eff <- se$d_hat; d_eff[is.na(d_eff)] <- 0
  
  # 缺失填补
  X <- impute_mode_matrix(X)
  
  # 逐位点贡献（不做任何等位基因对齐/翻转）
  contrib <- numeric(nrow(X))
  for (j in seq_len(ncol(X))) {
    g  <- X[, j]
    a  <- a_eff[j]
    d  <- d_eff[j]
    pd <- se$p_d[j]
    
    c_add <- numeric(length(g))
    c_add[g == 0L] <- -a
    c_add[g == 1L] <-  0
    c_add[g == 2L] <-  a
    
    if (!is.na(pd) && pd < dominance_p_threshold) {
      c_dom <- numeric(length(g)); c_dom[g == 1L] <- d
      contrib <- contrib + c_add + c_dom
    } else {
      contrib <- contrib + c_add
    }
  }
  prs_vec <- prs_vec + contrib
}

# 输出与评估
prs_a_d_noclump_path <- file.path(pheno_dir, "testing_prs_a_d_causal.txt")
fwrite(data.table(IID = test_ids, PRS = prs_vec), prs_a_d_noclump_path, sep = "\t")

ph_test <- fread(output_testing_pheno)[, .(IID, pheno)]
prs_vec <- fread(file.path(pheno_dir, "testing_prs_a_d_causal.txt"))
res_nc  <- merge(prs_vec, ph_test, by = "IID")
r_nc    <- suppressWarnings(cor(res_nc$pheno, res_nc$PRS, use = "complete.obs"))
cat(sprintf("[TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=%d  cor=%.6f\n", nrow(res_nc), r_nc))
cat("Saved PRS_a_d (no clump, causal-only, no allele alignment) -> ",
    prs_a_d_noclump_path, "\n", sep = "")
#[TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.706299
#[TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.703325
#[TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.695476
#100 comdom
#[TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.691878
#50comdom
# [TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.696582
# 100 add
# [TEST causal ALL-SNP a+d (no clump, no allele alignment)] N=20000  cor=0.702076