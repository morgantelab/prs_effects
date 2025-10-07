#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse command-line arguments
rep <- as.integer(args[1])
out_dir <- args[2]  # output directory
test_ids_file <- args[3]
effects_file <- args[4]
geno_pattern <- args[5]
plink_path <- args[6]

cat("\nStarting PRS calculation for replicate", rep, "\n")
cat("Output directory:", out_dir, "\n")
cat("Test ID file:", test_ids_file, "\n")
cat("Effect file:", effects_file, "\n")
cat("Geno pattern:", geno_pattern, "\n")
cat("PLINK path:", plink_path, "\n")

library(data.table)
library(genio)

# Create temp directory
temp_dir <- file.path(out_dir, sprintf("temp_rep%d", rep))
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created temp directory:", temp_dir, "\n")
}

# Load test set IDs
test_ids <- fread(test_ids_file, header = FALSE)[[1]]

# Load effect estimates
eff <- fread(effects_file)
if ("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")

# Initialize PRS vector
prs_vector <- numeric(length(test_ids))

# Safe addition function: 保留缺失但不影响其他 SNP
safe_add <- function(x, y) {
  out <- x
  both <- !is.na(x) & !is.na(y)
  only_y <-  is.na(x) & !is.na(y)
  out[both] <- x[both] + y[both]  # 两个都有值就加
  out[only_y] <- y[only_y]        # x 缺失但 y 有值 → 用 y
  out
}

# Process each chromosome
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  chr_snps <- eff[CHR == chr, .(SNP, ref_allele, alt_allele, a_hat, d_hat)]
  if (nrow(chr_snps) == 0) {
    cat("  No SNPs for chromosome", chr, "\n")
    next
  }
  
  snp_list_file <- file.path(temp_dir, sprintf("chr%d_snps.txt", chr))
  fwrite(chr_snps[, .(SNP)], snp_list_file, col.names = FALSE)
  
  plink_file <- paste0(geno_pattern, chr, "_filtered")
  temp_prefix <- file.path(temp_dir, sprintf("chr%d", chr))
  
  plink_cmd <- sprintf(
    "%s --bfile %s --keep %s --extract %s --make-bed --out %s --memory 2000 --threads 40",
    plink_path, plink_file, test_ids_file, snp_list_file, temp_prefix
  )
  system(plink_cmd)
  
  plink_data <- genio::read_plink(temp_prefix, verbose = FALSE)
  snp_ids <- plink_data$bim$id
  X <- t(plink_data$X)  # 个体为行，SNP 为列
  
  # 对齐 SNP 顺序
  snp_effects <- chr_snps[match(snp_ids, chr_snps$SNP), ]
  
  prs_contrib <- numeric(nrow(X))
  
  # 循环 SNP 计算贡献
  for (j in seq_along(snp_ids)) {
    g <- X[, j]   # 保留 NA
    a <- snp_effects$a_hat[j]
    d <- snp_effects$d_hat[j]
    
    contrib <- rep(NA_real_, length(g))
    contrib[g == 0] <- -a
    contrib[g == 1] <- d
    contrib[g == 2] <- a
    
    # 累加到总 PRS，使用 safe_add 保留 NA 但不污染其他
    prs_contrib <- safe_add(prs_contrib, contrib)
  }
  
  # 累加到最终 PRS
  prs_vector <- safe_add(prs_vector, prs_contrib)
}

unlink(temp_dir, recursive = TRUE)

# 保存结果
prs_df <- data.table(IID = test_ids, PRS = prs_vector)
output_file <- file.path(out_dir, sprintf("test_prs_scan_rep%d.txt", rep))
fwrite(prs_df, output_file, sep = "\t")
cat("Saved PRS results to:", output_file, "\n")