suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
})

## Files produced earlier are assumed to exist:
## - pheno_dir
## - plink2_path
## - geno_pattern (prefix "ukb_chr" -> "ukb_chr{CHR}_filtered")
## - output_training_ID, output_testing_ID
## - output_training_pheno, output_testing_pheno
## - file.path(pheno_dir, "clumped_snps.snplist")

clumped_snplist_file <- file.path(pheno_dir, "clumped_snps.snplist")
if (!file.exists(clumped_snplist_file)) stop("Missing clumped SNP list: ", clumped_snplist_file)

train_ids <- fread(output_training_ID, header = FALSE, col.names = c("FID","IID"))
test_ids  <- fread(output_testing_ID,  header = FALSE, col.names = c("FID","IID"))
training_pheno <- fread(
  output_training_pheno,
  colClasses = c(FID="character", IID="character", pheno="numeric")
)
# keep IDs as character to match FAM/IID
training_pheno[, IID := as.character(IID)]

test_pheno  <- fread(output_testing_pheno)

## --- helper: estimate a,d per chromosome on TRAINING set ---
calculate_welch <- function(data, snp_cols) {
  n_snps <- length(snp_cols)
  results <- data.frame(
    SNP   = snp_cols,
    a_hat = rep(NA_real_, n_snps),
    d_hat = rep(NA_real_, n_snps),
    p_a   = rep(NA_real_, n_snps),
    p_d   = rep(NA_real_, n_snps)
  )
  
  pheno <- data$pheno
  
  for (j in seq_len(n_snps)) {
    snp <- snp_cols[j]
    X <- data[[snp]]
    
    # 基因型分组
    AA <- X == 2
    Aa <- X == 1
    aa <- X == 0
    
    n_AA <- sum(AA, na.rm = TRUE)
    n_Aa <- sum(Aa, na.rm = TRUE)
    n_aa <- sum(aa, na.rm = TRUE)
    
    
    # 组内表型
    y_AA <- pheno[AA]
    y_Aa <- pheno[Aa]
    y_aa <- pheno[aa]
    
    mean_AA <- mean(y_AA, na.rm = TRUE)
    mean_Aa <- mean(y_Aa, na.rm = TRUE)
    mean_aa <- mean(y_aa, na.rm = TRUE)
    var_AA  <- var(y_AA,  na.rm = TRUE)
    var_Aa  <- var(y_Aa,  na.rm = TRUE)
    var_aa  <- var(y_aa,  na.rm = TRUE)
    
    # Additive effect (Welch)
    a_hat <- (mean_AA - mean_aa) / 2
    var_a <- 0.25 * (var_AA / n_AA + var_aa / n_aa)
    if (var_a <= 0) next
    # Welch-Satterthwaite df
    df_a <- (var_AA/n_AA + var_aa/n_aa)^2 /
      ((var_AA/n_AA)^2/(n_AA-1) + (var_aa/n_aa)^2/(n_aa-1))
    t_a <- a_hat / sqrt(var_a)
    p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)
    
    # Dominance effect (Welch)
    d_hat <- mean_Aa - (mean_AA + mean_aa)/2
    var_d <- var_Aa/n_Aa + 0.25*(var_AA/n_AA + var_aa/n_aa)
    
    df_d <- (var_d)^2 / (
      (var_Aa/n_Aa)^2/(n_Aa-1) +
        (var_AA/(4*n_AA))^2/(n_AA-1) +
        (var_aa/(4*n_aa))^2/(n_aa-1)
    )
    t_d <- d_hat / sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}

## ========= 1) SNPs to evaluate: GWAS-clumped/thresholded =========
snps_keep <-  file.path(pheno_dir, "clumped_snps.snplist")

chunk_size <- 5000 

all_results <- list()
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  # 输入前缀（你的已过滤 bfile 前缀）
  full_plink <- paste0(geno_pattern, chr, "_filtered")
  if (!file.exists(paste0(full_plink, ".bim"))) {
    cat("Warning: No .bim file for chromosome", chr, "- skipping\n")
    next
  }
  
  # 临时前缀（只保留训练集个体）
  temp_prefix <- file.path(tempdir(), paste0("chr", chr, "_train_all"))
  temp_bed <- paste0(temp_prefix, ".bed")
  temp_bim <- paste0(temp_prefix, ".bim")
  temp_fam <- paste0(temp_prefix, ".fam")
  
  # 仅保留训练集样本以及通过GWAS的C+T的SNP
  plink_cmd <- sprintf(
    "%s --bfile %s --keep %s --keep-allele-order --extract %s --make-bed --out %s --threads %d",
    shQuote(plink2_path),
    shQuote(full_plink),
    shQuote(output_training_ID),
    shQuote(snps_keep),
    shQuote(temp_prefix),
    20
  )
  cat("Running:", plink_cmd, "\n")
  system(plink_cmd)
  
  if (!file.exists(temp_bed)) {
    cat("Warning: Plink2 failed for chromosome", chr, "- skipping\n")
    next
  }
  
  # FAM（训练集样本顺序）
  fam <- fread(temp_fam, col.names = c("FID","IID","pat","mat","sex","pheno"))
  
  # BIM（位点注释）
  bim_tmp <- fread(
    temp_bim,
    col.names = c("chr","SNP","pos_cM","pos","alt","ref")
  )
  
  # 按 FAM 顺序对齐训练表型
  pheno_vector <- training_pheno$pheno[match(fam$IID, training_pheno$IID)]
  
  # 训练集基因型：个体 × 全部位点
  bed <- BEDMatrix(temp_bed)  # 列名为 SNP ID（来自 BIM）
  
  n_snps_total <- ncol(bed)
  if (is.null(n_snps_total) || n_snps_total == 0) {
    cat("No genotype data loaded for chr", chr, "\n")
    unlink(c(
      paste0(temp_prefix, ".bed"),
      paste0(temp_prefix, ".bim"),
      paste0(temp_prefix, ".fam"),
      paste0(temp_prefix, ".log"),
      paste0(temp_prefix, ".nosex")
    ), force = TRUE)
    next
  }
  
  chr_results <- vector("list", ceiling(n_snps_total / chunk_size))
  idx <- 1L
  
  for (start_col in seq(1L, n_snps_total, by = chunk_size)) {
    end_col <- min(start_col + chunk_size - 1L, n_snps_total)
    snp_idx  <- start_col:end_col
    
    # ↓↓↓ 用 BIM 里的纯 SNP ID，而不是 BEDMatrix 的列名（可能带等位基因）
    snp_ids_block <- bim_tmp$SNP[snp_idx]
    
    # 提取该块基因型矩阵（训练集个体）
    snp_block <- bed[, snp_idx, drop = FALSE]
    
    # ↓↓↓ 把列名改成纯 SNP ID，保证与 BIM 一致
    colnames(snp_block) <- snp_ids_block
    
    # 组成数据框：pheno + SNP 列（保持原 SNP 名称）
    df_block <- data.frame(
      pheno = pheno_vector,
      as.data.frame(snp_block),
      check.names = FALSE
    )
    
    # ↓↓↓ 传入的 snp 列表也用纯 SNP ID
    res_block <- calculate_welch(df_block, snp_ids_block)
    
    chr_results[[idx]] <- res_block
    cat(sprintf("  Chr %d: processed SNPs %d-%d / %d\n",
                chr, start_col, end_col, n_snps_total))
    idx <- idx + 1L
  }
  
  # 合并该染色体结果并加注释
  results_dt <- rbindlist(chr_results, use.names = TRUE, fill = TRUE)
  
  results_final <- merge(
    results_dt,
    bim_tmp[, .(SNP, pos, ref, alt)],
    by = "SNP",
    all.x = TRUE
  )[, .(
    chromosome = chr,
    position   = pos,
    SNP,
    ref_allele = ref,
    alt_allele = alt,
    a_hat, d_hat, p_a, p_d
  )]
  
  all_results[[length(all_results) + 1L]] <- results_final
  
  # 清理临时文件
  unlink(c(
    paste0(temp_prefix, ".bed"),
    paste0(temp_prefix, ".bim"),
    paste0(temp_prefix, ".fam"),
    paste0(temp_prefix, ".log"),
    paste0(temp_prefix, ".nosex")
  ), force = TRUE)
  
  cat("Finished chromosome", chr, "with", nrow(results_final), "SNPs.\n")
}

# 写出训练集效应表
effects_out <- file.path(pheno_dir, "effects_training_method2.txt")
fwrite(rbindlist(all_results, use.names = TRUE, fill = TRUE),
       effects_out, sep = "\t")
cat("✅ Wrote effects → ", effects_out, "\n", sep = "")


## --- Calculate PRS_a_d in testing sets ---
library(genio)
stopifnot(file.exists(effects_out))

# Load Welch effects estimated on TRAINING set
effects <- fread(effects_out)
stopifnot(all(c("SNP","a_hat","d_hat","p_d") %in% names(effects)))
effects[, chromosome := as.integer(chromosome)]
# Testing IDs & pheno
test_ids  <- fread(output_testing_ID,  header = FALSE, col.names = c("FID","IID"))
test_pheno <- fread(output_testing_pheno)
test_pheno[, IID := as.character(IID)]

# Testing IDs (assume two-column keep: FID IID)
test_ids <- fread(output_testing_ID, header = FALSE)[[2]]
tmp_dir <- file.path(pheno_dir, "tmp_prs_ad"); dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
dominance_p_threshold <- 0.01
# 按 SNP 众数填补缺失 (0/1/2)
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
prs_vec <- numeric(length(test_ids)); names(prs_vec) <- test_ids

for (chr in 1:22) {
  se_chr <- effects[chromosome == chr, .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_d)]
  if (!nrow(se_chr)) next
  
  # 提取该染色体测试集基因型（仅 clumped SNP）
  snplist_chr <- file.path(tmp_dir, sprintf("prs_chr%d.snps", chr))
  fwrite(se_chr[, .(SNP)], snplist_chr, col.names = FALSE)
  
  bfile      <- sprintf("%s%d_filtered", geno_pattern, chr)
  out_prefix <- file.path(tmp_dir, sprintf("prs_extract_chr%d", chr))
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
  X       <- t(gd$X)  # N_ind x N_snp
  
  # 行顺序对齐
  idx <- match(test_ids, fam_ids)
  X   <- X[idx, , drop = FALSE]
  
  # 列顺序对齐（SNP 顺序）
  se <- se_chr[match(snp_ids, se_chr$SNP)]
  keep <- !is.na(se$SNP)
  X  <- X[, keep, drop = FALSE]
  se <- se[keep]
  snp_ids <- snp_ids[keep]
  
  a_eff <- se$a_hat
  d_eff <- se$d_hat
  
  # 缺失用众数填补
  X <- impute_mode_matrix(X)
  
  # 逐 SNP 累加：g=0→-a, g=1→(0 或 +d), g=2→+a
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

prs_a_d_path <- file.path(pheno_dir, "testing_prs_a_d_method2.txt")
fwrite(data.table(IID = test_ids, PRS = prs_vec), prs_a_d_path, sep = "\t")
ph_test <- fread(output_testing_pheno)[, .(IID, pheno)]
prs_vec <- fread(prs_a_d_path)
res     <- merge(prs_vec, ph_test, by = "IID")
rval    <- suppressWarnings(cor(res$pheno, res$PRS, use = "complete.obs"))
cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
cat("Saved PRS_a_d -> ", prs_a_d_path, "\n", sep = "")
#100 add
# > cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
# [TEST clumped a+d] N=20000  cor=0.539659
# 100 comdom
# > cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
# [TEST clumped a+d] N=20000  cor=0.517488
# 50com50add
# > cat(sprintf("[TEST clumped a+d] N=%d  cor=%.6f\n", nrow(res), rval))
# [TEST clumped a+d] N=20000  cor=0.531608
