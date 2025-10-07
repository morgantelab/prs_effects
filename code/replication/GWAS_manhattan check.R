library(data.table)
library(ggplot2)
library(ggrepel)

results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2"
gwas_dir <- file.path(results_dir, "gwas_test")  # New directory for GWAS results
causal_snps_file <- file.path(results_dir, "pheno/causal_snps1.txt")
causal_snps <- scan(causal_snps_file, what = "character", quiet = TRUE)
# 构建文件路径
rep = 1
gwas_files <- file.path(gwas_dir, sprintf("gwas_rep%d_chr%d.pheno.glm.linear", rep, 1:22))

# 读取并合并
gwas_list <- lapply(1:22, function(chr) {
  f <- gwas_files[chr]

  dat <- fread(f)  # 跳过第一行的 ## 注释
  setnames(dat, toupper(names(dat)))  # 确保列名大写
  
  # 提取需要的列
  dat[, `:=`(CHROM = as.numeric(`#CHROM`), 
             POS = as.numeric(POS), 
             P = as.numeric(P), 
             BETA = as.numeric(BETA), 
             ID = as.character(ID))]
  
  return(dat[!is.na(P), .(ID, CHROM, POS, P, BETA)])
})

# 合并所有染色体
gwas_all <- rbindlist(gwas_list)
setkey(gwas_all, ID)

gwas_all[, is_causal := ID %in% causal_snps]

manhattan <- gwas_all[order(CHROM, POS)]
manhattan[, cumpos := c(0, cumsum(head(diff(range(POS)), -1)))[CHROM], by = CHROM]
manhattan[, x := cumpos + POS]

# 绘图
ggplot(manhattan, aes(x = x, y = -log10(P))) +
  # Points for non-causal SNPs (alternating black/gray by chromosome)
  geom_point(
    data = subset(manhattan, !is_causal),
    aes(color = factor(CHROM %% 2)),  # Alternate colors by even/odd chr
    size = 0.8,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("black", "gray50")) +  # Black (odd chr) & Gray (even chr)
  
  # Highlight causal SNPs in red
  geom_point(
    data = subset(manhattan, is_causal),
    color = "red",
    size = 1.5
  ) +
  
  # Label causal SNPs
  ggrepel::geom_text_repel(
    data = subset(manhattan, is_causal),
    aes(label = ID),
    size = 3,
    nudge_y = 1,
    point.padding = 0.3,
    segment.alpha = 0.5,
    max.overlaps = 20
  ) +
  
  # Genome-wide significance line
  geom_hline(
    yintercept = -log10(5e-8),
    linetype = "dashed",
    color = "blue",
    alpha = 0.8
  ) +
  
  # Theme and axis formatting
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Chromosome") +
  ylab("-log10(P)") +
  scale_x_continuous(
    breaks = sapply(1:22, function(chr) {
      mean(subset(manhattan, CHROM == chr)$x)
    }),
    labels = 1:22
  ) +
  ggtitle(paste0("Manhattan Plot (rep = ", rep, ") - Causal SNPs Labeled"))


# 重新生成 true a 向量（确保顺序一致！）
set.seed(1001)
a <- rnorm(1000, mean = 0, sd = sqrt(2))

# 生成 dom_type 和 d（虽然这里只用 a）
dom_type <- rep(c("complete", "incomplete", "additive"), times = c(0, 200, 800))
d <- numeric(1000)
additive_idx <- which(dom_type == "additive")
d[additive_idx] <- 0

# 关键：必须确保 `a` 的顺序与 `causal_snps` 文件中 SNP 的顺序一致
# 假设 causal_snps 是按生成顺序保存的，则：
names(a) <- causal_snps  # 把 rsID 作为名字

# 提取因果 SNP 的 GWAS beta
causal_results <- gwas_all[ID %in% causal_snps, .(ID, BETA)]
causal_results[, true_a := a[match(ID, names(a))]]

# 检查是否全部匹配
if (any(is.na(causal_results$true_a))) {
  warning("Some causal SNPs not found in true effect vector!")
}

# 绘图
ggplot(causal_results, aes(x = true_a, y = BETA)) +
  geom_point(alpha = 0.8, color = "blue", size = 2.5) +
  xlab("True Additive Effect (a)") +
  ylab("Estimated GWAS Beta") +
  theme_bw() +
  ggtitle("GWAS Beta vs True Effect for Causal SNPs") 
