#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# ---------------------- 命令行参数 ----------------------
rep           <- as.integer(args[1])
geno_pattern  <- args[2]
chunk_size    <- as.integer(args[3])
output_effects<- args[4]
output_train  <- args[5]
output_test   <- args[6]
pheno_file    <- args[7]
seed_val      <- as.integer(args[8])
num_train     <- as.integer(args[9])
num_sample    <- as.integer(args[10])
plink2_path   <- args[11]   # 注意：已移到第 11 个参数（删除了 causal_snps 参数）

cat("\nConverted values:\n")
cat("Geno file pattern:", geno_pattern, "\n")
cat("rep:", rep, "\n")
cat("seed_val:", seed_val, "\n")
cat("chunk_size:", chunk_size, "\n")
cat("output_effects:", output_effects, "\n")
cat("output_train:", output_train, "\n")
cat("output_test:", output_test, "\n")
cat("pheno_file:", pheno_file, "\n")
cat("Number of Training set:", num_train, "\n")
cat("Number of Total Sample size:", num_sample, "\n")
cat("plink2_path:", plink2_path, "\n")

# ---------------------- 依赖包 ----------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(genio)
  library(BEDMatrix)
})

# ---------------------- 函数：Welch 计算 ----------------------
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

# ---------------------- 读入表型并划分训练/测试 ----------------------
pheno <- fread(
  pheno_file,
  colClasses = c("FID" = "character", "IID" = "character", "pheno" = "numeric")
)

cat("Successfully read phenotype files for replicate", rep, "\n")

set.seed(seed_val)
sample_ids <- pheno$IID
if (length(sample_ids) < num_train) {
  stop("num_train is larger than available IDs in phenotype file.")
}
train_ids <- sample(sample_ids, num_train)
test_ids  <- setdiff(sample_ids, train_ids)

train_pheno <- pheno[IID %in% train_ids]
test_pheno  <- pheno[IID %in% test_ids]

# 写训练/测试 ID 文件（供 plink 的 --keep 使用）
fwrite(train_pheno[, .(FID, IID)], output_train, col.names = FALSE, sep = "\t", quote = FALSE)
fwrite(test_pheno[, .(FID, IID)],  output_test,  col.names = FALSE, sep = "\t", quote = FALSE)

cat("Created training/test ID files:\n",
    "- Training:", output_train, "\n",
    "- Test:", output_test, "\n")

# ---------------------- 按染色体处理（训练集个体 + 全部位点） ----------------------
all_results <- list()

for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  # 输入前缀
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
  
  # 运行 plink：保留训练集，保留所有 SNP
  plink_cmd <- sprintf(
    "%s --bfile %s --keep %s --keep-allele-order --make-bed --out %s --threads %d",
    plink2_path,
    full_plink,
    output_train,
    temp_prefix,
    8
  )
  cat("Running:", plink_cmd, "\n")
  system(plink_cmd, intern = FALSE)
  
  if (!file.exists(temp_bed)) {
    cat("Warning: Plink2 failed for chromosome", chr, "- skipping\n")
    next
  }
  
  # 读取 FAM（训练集样本顺序）
  fam <- fread(temp_fam, col.names = c("FID","IID","pat","mat","sex","pheno"))
  # 读取 BIM（所有位点信息）
  bim_tmp <- fread(temp_bim,
                   col.names = c("chr","SNP","pos_cM","pos","alt","ref"))
  
  # 训练集表型向量，按 fam 顺序对齐
  pheno_vector <- train_pheno$pheno[match(fam$IID, train_pheno$IID)]
  
  # 读取 BED（训练集个体 × 全部位点）
  bed <- BEDMatrix(temp_bed, simple_names = TRUE)
  
  # 分块处理列（SNP）
  n_snps_total <- ncol(bed)
  if (is.null(n_snps_total) || n_snps_total == 0) {
    cat("No genotype data loaded for chr", chr, "\n")
    # 清理临时文件
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
    snp_idx <- start_col:end_col
    snp_cols <- colnames(bed)[snp_idx]
    
    # 提取该块的基因型矩阵（已是训练集个体）
    snp_block <- bed[, snp_idx, drop = FALSE]
    
    # 构建数据框：pheno + SNP 列
    df_block <- data.frame(
      pheno = pheno_vector,
      as.data.frame(snp_block),
      check.names = FALSE
    )
    
    # 计算该块的效应
    res_block <- calculate_welch(df_block, snp_cols)
    chr_results[[idx]] <- res_block
    cat(sprintf("  Chr %d: processed SNPs %d-%d / %d\n",
                chr, start_col, end_col, n_snps_total))
    idx <- idx + 1L
  }
  
  # 合并当前染色体的块结果
  results_dt <- rbindlist(chr_results, use.names = TRUE, fill = TRUE)
  
  # 合并位点注释
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

# ---------------------- 汇总与写出 ----------------------
if (length(all_results) > 0) {
  final_dt <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  setorder(final_dt, chromosome, position)
  fwrite(final_dt, output_effects, sep = "\t")
  cat("✅ Final effects written to:", output_effects, "\n")
} else {
  cat("⚠️ No results to write for replicate", rep, "\n")
}

