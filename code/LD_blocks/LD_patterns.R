library(data.table)
library(ggplot2)
library(bigsnpr)

# Set paths
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/LD_decay"
geno_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/filter"
# Load LD blocks
ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr5.bed")

# Read BIM file
bim <- fread(paste0(geno_dir, "/ukb_chr5_filtered.bim"),
             select = c(2, 4), col.names = c("snp", "pos"))

# Process each block
all_ld_data <- list()

for (i in 1:nrow(ldetect_chr)) {
  start <- ldetect_chr$start[i]
  stop <- ldetect_chr$stop[i]

  # Construct PLINK command
  output_prefix <- file.path(results_dir, paste0("block_", i))
  system(paste(
    plink2_path,
    "--bfile", paste0(geno_dir, "/ukb_chr5_filtered"),
    "--chr 5",
    "--from-bp", start,
    "--to-bp", stop,
    "--r2-unphased",
    "--ld-window-r2 0",
    "--ld-window 100000",   # Large window to cover all possible pairs
    "--ld-window-kb 1000", # Large distance (100Mb) to cover chromosome
    "--out", output_prefix,
    "--memory 8000",  # Less memory per block
    "--threads 12"     # Fewer threads per block
  ))

  # 处理.vcor文件 - 直接读取配对结果
  vcor_file <- paste0(output_prefix, ".vcor")
  if (file.exists(vcor_file)) {
    # 读取文件并处理
    system(paste("cat", vcor_file, "| sed 1,1d | awk -F \"\t\" 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS=\"\t\"}{print abs($5-$2),$7}' >",
                 paste0(output_prefix, ".ld.summary")))

    # 读取摘要文件
    ld_summary <- fread(paste0(output_prefix, ".ld.summary"), header = FALSE, sep = "\t")
    colnames(ld_summary) <- c("dist", "rsq")

    # 添加区块标识
    ld_summary$block <- i

    # 保存到列表
    all_ld_data[[i]] <- ld_summary

    cat("  找到", nrow(ld_summary), "个SNP对\n")

    file.remove(vcor_file)
    file.remove(paste0(output_prefix, ".log"))
  }
}

# 合并所有结果
ld_combined <- rbindlist(all_ld_data)

# 保存原始数据
fwrite(ld_combined, file.path(results_dir, "chr5_all_ld_pairs.csv"))

# 创建摘要文件
cat("创建摘要文件...\n")
ld_summary_file <- file.path(results_dir, "chr5.ld.summary")
fwrite(ld_combined[, .(dist,  rsq)], ld_summary_file,
       sep = " ", col.names = FALSE, quote = FALSE)

# 保存结果
fwrite(ld_combined, file.path(results_dir, "chr5_ld_results.csv"))


###--------------One whole chromosome-----------------

library(data.table)
library(ggplot2)

# Set paths
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/LD_decay"
geno_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/filter"

# Output prefix for whole chromosome
output_prefix <- file.path(results_dir, "whole_chr1")

# Run PLINK for entire chromosome 1
system(paste(
  plink2_path,
  "--bfile", paste0(geno_dir, "/ukb_chr1_filtered"),
  "--chr 1",
  "--r2-unphased",
  "--ld-window-r2 0",
  "--ld-window 100000",   # Large window to cover all possible pairs
  "--ld-window-kb 1000", # Large distance (100Mb) to cover chromosome
  "--out", output_prefix,
  "--memory 8000",
  "--threads 24"
))

# Process .vcor file to calculate distances and r²
vcor_file <- paste0(output_prefix, ".vcor")

vcor_data <- fread(vcor_file, sep = "\t", header = TRUE)


# Read and process results
ld_data <- vcor_data[, .(
  dist = abs(POS_B - POS_A),  # 计算绝对距离
  rsq = UNPHASED_R2                  # 直接提取R2列
)]


# # Save full results
# fwrite(ld_data, file.path(results_dir, "chr1_full_ld_results.csv"))
# 
# # Cleanup temporary files
# file.remove(paste0(output_prefix, ".log"))
# 
# ld_data <-fread(file.path(results_dir, "chr1_full_ld_results.csv"))

# Create binned summary for plotting
bin_size <- 10000  # 10kb bins
ld_data[, bin := floor(dist / bin_size) * bin_size]
# 创建包含平均值和中位数的摘要
ld_summary_comparison <- ld_data[, .(
  mean_rsq = mean(rsq, na.rm = TRUE),
  median_rsq = median(rsq, na.rm = TRUE),
  q25 = quantile(rsq, 0.25, na.rm = TRUE),
  q75 = quantile(rsq, 0.75, na.rm = TRUE)
), by = bin][order(bin)]

# 绘制对比图
ggplot(ld_summary_comparison, aes(x = bin)) +
  geom_line(aes(y = mean_rsq, color = "Mean"), linewidth = 1) +
  geom_line(aes(y = median_rsq, color = "Median"), linewidth = 1) +
  scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) +
  labs(x = "Distance (bp)", y = "r²", 
       title = "LD Decay Comparison",
       color = "Statistic") +
  scale_y_continuous(limits = c(0, max(ld_summary_comparison$mean_rsq)), expand = c(0, 0)) +
  scale_x_continuous(breaks=c(0,1E5,2E5,4E5,6E5),
                     labels=c("0","1","2","4","6")) +
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.9),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


ggsave(file.path(results_dir, "chr1_ld_decay_plot.png"), 
       width = 10, height = 6, dpi = 300)

### 在unfiltered数据上验证####

# Set paths
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/LD_decay"
bed_bim_dir <- "/data2/morgante_lab/data/ukbiobank/genotypes/array"
fam_dir <- "/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024"

# Output prefix
output_prefix <- file.path(results_dir, "unfiltered_chr1_ld")

# Chromosome to process
chr <- 1

# Construct file paths
bed_file <- file.path(bed_bim_dir, paste0("ukb22418_c", chr, "_b0_v2_s488243.bed"))
bim_file <- file.path(bed_bim_dir, paste0("ukb22418_c", chr, "_b0_v2_s488243.bim"))
fam_file <- file.path(fam_dir, paste0("ukb22418_c", chr, "_b0_v2_s487945.fam"))

# Define distance bins (up to 1Mb)
bins <- seq(0, 1000000, by = 10000)  # 10kb bins
bin_centers <- bins[-1] - 5000  # 计算每个bin的中心位置

# Initialize results storage
block_results <- list()

# Process each LD block
for (i in 1:nrow(ldetect_chr)) {
  start <- ldetect_chr$start[i]
  stop <- ldetect_chr$stop[i]
  
  # Output prefix for this block
  output_prefix <- file.path(results_dir, paste0("block_", i))
  
  # Build PLINK command for this block
  plink_cmd <- paste(
    plink2_path,
    "--bed", bed_file,
    "--bim", bim_file,
    "--fam", fam_file,
    "--rm-dup exclude-all",
    "--snps-only just-acgt",
    "--maf 0.01",
    "--mac 5",
    "--geno 0.1",
    "--hwe 1e-5",
    "--chr", chr,
    "--from-bp", start,
    "--to-bp", stop,
    "--r2-unphased",
    "--ld-window-r2 0",
    "--ld-window 100000",
    "--ld-window-kb 1000",
    "--out", output_prefix,
    "--memory 8000",
    "--threads 12"
  )
  
  # Run PLINK for this block
  system(plink_cmd)
  
  # Process output
  vcor_file <- paste0(output_prefix, ".vcor")
  
  
  # Read .vcor file
  vcor_data <- fread(vcor_file, sep = "\t", header = TRUE)
  
  # Calculate distances
  vcor_data[, dist := abs(POS_B - POS_A)]
  # Bin distances and calculate mean r² per bin
  vcor_data[, bin_index := findInterval(dist, bins)]
  # Calculate mean r² for each bin
  bin_summary <- vcor_data[, .(mean_rsq = mean(UNPHASED_R2, na.rm = TRUE)), 
                           by = bin_index]
  # Add block information
  bin_summary[, block := i]
  bin_summary[, start_bp := start]
  bin_summary[, end_bp := stop]
  
  # Save to results
  block_results[[i]] <- bin_summary
  
  # Clean up
  file.remove(vcor_file)
  # Progress update
  if (i %% 10 == 0) message("Processed ", i, " of ", nrow(ldetect_chr), " blocks")
}

# Combine all block results
all_blocks <- rbindlist(block_results)

# Map bin indices to bin centers
all_blocks[, bin_center := bins[bin_index] + 5000]

# Generate LD decay plot with one line per block
ggplot() +
  # Individual block lines (faint)
  geom_line(data = all_blocks, 
            aes(x = bin_center, y = mean_rsq, group = block),
            alpha = 0.1, color = "blue") +
  scale_x_continuous(
    name = "Distance (bp)",
    breaks = seq(0, 1000000, by = 100000),
    labels = scales::comma
  ) +
  scale_y_continuous(
    name = "Mean r²",
    limits = c(0, 0.3)
  ) +
  labs(title = "LD Decay by Block on Chromosome 1") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


###---------------Summary stats-------------------------------###
# 加载必要包
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# 设置路径 (根据实际修改)
ldetect_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr1.bed"
bim_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/filter/ukb_chr1_filtered.bim"


# 读取并预处理数据
ld_blocks <- fread(ldetect_file, col.names = c("chr", "start", "end")) %>%
  mutate(chr = as.character(chr),
         block_id = paste0("block_", 1:n()),
         start_1based = start + 1)  # 转换为1-based坐标

snp_positions <- fread(bim_file, 
                       col.names = c("chr", "snp_id", "cM", "pos", "allele1", "allele2")) %>%
  select(chr, snp_id, pos) %>%
  mutate(chr = as.character(chr))


# 修复方法2：更高效的非rowwise方法
snp_assignments <- ld_blocks %>%
  mutate(block_chr = chr) %>%  # 创建副本避免冲突
  inner_join(snp_positions, by = character(0)) %>%  # 创建所有可能的组合
  filter(pos >= start_1based,
         pos < end) %>%
  select(block_id, block_start = start, block_end = end, 
         snp_id, snp_pos = pos)

# 计算每个区块内SNP间距统计量（使用修复后的数据）
block_distances <- snp_assignments %>%
  group_by(block_id, block_start, block_end) %>%
  arrange(snp_pos) %>%
  mutate(distance = c(NA, diff(snp_pos))) %>%
  filter(!is.na(distance)) %>%
  summarise(
    n_snps = n() + 1,  # 包括第一个SNP
    min_dist = min(distance),
    max_dist = max(distance),
    mean_dist = mean(distance),
    median_dist = median(distance),
    q25 = quantile(distance, 0.25),
    q75 = quantile(distance, 0.75),
    .groups = "drop"
  ) %>%
  # 添加区块长度
  mutate(block_length = block_end - block_start)


write.table(block_distances, paste0(results_dir, "CHR1_LD_blocks_summstats.xlsx"),
            quote = FALSE, row.names = F, col.names = T)

# 全局统计
global_stats <- block_distances %>%
  summarise(
    total_blocks = n(),
    blocks_with_snps = sum(n_snps > 1),
    mean_snps_per_block = mean(n_snps),
    mean_distance = mean(mean_dist),
    median_distance = median(median_dist),
    min_distance = min(min_dist),
    max_distance = max(max_dist)
  )


# Save plots
ggsave("chr1_distance_distribution.png", plot1, width = 8, height = 6, dpi = 300)
ggsave("chr1_length_vs_distance.png", plot2, width = 8, height = 6, dpi = 300)
                      
# 输出结果
cat("=== 染色体1 LD区块统计 ===\n")
cat(sprintf("总区块数: %d\n", global_stats$total_blocks))
cat(sprintf("包含SNP的区块: %d (%.1f%%)\n", 
            global_stats$blocks_with_snps,
            global_stats$blocks_with_snps/global_stats$total_blocks*100))
cat(sprintf("平均每个区块SNP数: %.1f\n", global_stats$mean_snps_per_block))
cat(sprintf("SNP间平均距离: %.0f bp\n", global_stats$mean_distance))
cat(sprintf("SNP间距离中位数: %.0f bp\n", global_stats$median_distance))

# 保存结果
fwrite(block_distances, "chr1_ld_block_distances.csv")
fwrite(global_stats, "chr1_global_ld_stats.csv")

# 可视化
# Plot 1: Distribution of mean SNP distances per LD block
ggplot(block_distances, aes(x = mean_dist)) +
  geom_histogram(fill = "orange", bins = 50, alpha = 0.8) +
  labs(title = "Distribution of Mean SNP Distances in LD Blocks (Chr1)",
       x = "Mean Distance Between SNPs (bp)", 
       y = "Number of LD Blocks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Plot 2: Relationship between block length and SNP distances
ggplot(block_distances, aes(x = block_length, y = mean_dist)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_log10(labels = scales::comma) +  # Format axis labels
  scale_y_log10(labels = scales::comma) +
  labs(title = "Relationship Between LD Block Length and SNP Distance",
       x = "LD Block Length (bp, log10 scale)", 
       y = "Mean SNP Distance (bp, log10 scale)",
       caption = "Red line shows linear regression fit") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 1, face = "italic"))

