library(data.table)
library(ggplot2)
library(bigsnpr)
library(dplyr)
library(stringr)

results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/LD_decay"

ld_combined <- fread(paste0(results_dir, "/chr1_ld_results.csv"))
# 生成LD衰减图 - 所有区块在同一图中
# 创建绘图
# 
ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr1.bed")
# 
# 
# ld_plot <- ggplot(ld_combined, aes(x = distance, y = R2, group = block)) +
#   geom_smooth(aes(color = factor(block)), method = "loess", 
#               se = FALSE, size = 0.3, alpha = 0.5) +
#   geom_smooth(method = "loess", color = "black", size = 1.5, se = FALSE) +  # 总体平均线
#   scale_x_log10(
#     breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
#     labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb")
#   ) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(
#     x = "Distance (bp)",
#     y = expression(LD~(r^2)),
#     title = "LD patterns in Chromosome 1",
#     subtitle = paste("Total of", nrow(ldetect_chr), "LD blocks |", nrow(ld_combined), "SNP pairs")
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 12),
#     legend.position = "none"
#   ) +
#   guides(color = "none")
# 
# 
# # 保存图像
# ggsave(file.path(results_dir, "chr1_ld_decay_all_blocks.pdf"), plot =  ld_plot, device = "pdf")
# 
# 
# selected_block <- 12
# ld_single_block <- ld_combined[ld_combined$block == selected_block, ]
# 
# ld_plot_single <- ggplot(ld_single_block, aes(x = distance, y = R2)) +
#   geom_smooth(method = "loess", color = "blue", size = 1, se = FALSE) +
#   scale_x_log10(
#     breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
#     labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb")
#   ) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(
#     x = "Distance (bp)",
#     y = expression(LD~(r^2)),
#     title = paste("LD Decay in Block", selected_block, "of Chromosome 1"),
#     subtitle = paste(nrow(ld_single_block), "SNP pairs")
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid.minor = element_blank(),
#     text = element_text(size = 12)
#   )
# 
# ggsave(file.path(results_dir, paste0("chr1_ld_decay_block_", selected_block, ".pdf")), 
#        plot = ld_plot_single, 
#        device = "pdf")

# 创建距离分箱
max_dist <- max(ld_combined$dist, na.rm = TRUE)
bin_size <- 10000 
breaks <- seq(from = 0, to = max_dist + bin_size, by = bin_size)

# 为所有数据添加分箱
ld_combined$dist_bin <- cut(ld_combined$dist, breaks = breaks, include.lowest = TRUE)

# 计算每个区块每个分箱的平均R2
block_binned <- ld_combined %>% 
  group_by(block, dist_bin) %>% 
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  ) 

# 计算分箱中点（使用 dist_bin 列）
block_binned$mid <- as.numeric(block_binned$dist_bin) * bin_size - bin_size/2

# 创建包含所有区块的LD衰减图
ld_plot <- ggplot() +
  # 绘制每个区块的线（半透明）
  geom_line(
    data = block_binned,
    aes(x = mid, y = mean_R2, group = factor(block)),
    alpha = 0.3, 
    color = "steelblue",
    size = 0.5
  ) +
  # 添加坐标轴标签
  labs(
    x = "Distance (bp)",
    y = expression(LD~(r^{2})),
    title = "LD Decay on Chromosome 1 by Block",
    subtitle = paste(n_distinct(block_binned$block), "LD blocks |", 
                     format(nrow(ld_combined), big.mark = ","), "SNP pairs")
  ) +
  scale_x_continuous(breaks=c(0,1E5,2E5,4E5,6E5,8E5),
                     labels=c("0","100kb","200kb","400kb","600kb","800kb"))+
  
  # 设置Y轴范围
  scale_y_continuous(
    limits = c(0, max(block_binned$mean_R2, na.rm = TRUE) * 1.1)
  ) +
  # 使用简洁主题
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),   # 去除次网格线
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )
ld_plot


ggsave(file.path(results_dir, "chr1_ld_decay_all_blocks.png"), plot =  ld_plot, device = "png")




### Whole chromosome
# # Set bin size to 10kb for balanced resolution
# bin_size <- 100  # 10kb bins
# 
# # Calculate maximum distance and create bins
# max_dist <- max(ld_combined$dist, na.rm = TRUE)
# breaks <- seq(0, max_dist + bin_size, by = bin_size)
# ld_combined$dist_bin <- cut(ld_combined$dist, breaks = breaks, include.lowest = TRUE)
# 
# # Compute binned averages
# whole_chr_binned <- ld_combined %>%
#   group_by(dist_bin) %>%
#   summarise(
#     mean_R2 = mean(R2, na.rm = TRUE),
#     n_pairs = n(),
#     .groups = "drop"
#   ) %>%
#   filter(!is.na(dist_bin)) %>%  # Remove NA bins
#   mutate(mid = (as.numeric(dist_bin) - 0.5) * bin_size)  # Bin midpoint
# 
# # Create LD decay plot for entire chromosome
# ld_plot_whole <- ggplot(whole_chr_binned, aes(x = mid, y = mean_R2)) +
#   geom_line(color = "red", size = 0.8) +  # Solid line for LD decay
#   scale_x_log10(
#     breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
#     labels = c("1kb", "10kb", "100kb", "1Mb", "10Mb"),
#     limits = c(1000, max_dist)  # Focus on relevant distances
#   ) +
#   scale_y_continuous(
#     limits = c(0, 1),  # Adjusted based on typical LD values
#     expand = c(0, 0)
#   ) +
#   labs(
#     x = "Distance (bp)",
#     y = expression(LD~(r^{2})),
#     title = "LD Decay on Chromosome 1 "
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12, face = "bold"),
#     plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5)
#   )
# 
# # Save high-resolution plot
# ggsave(
#   file.path(results_dir, "chr1_ld_decay_whole_chromosome.png"), 
#   plot = ld_plot_whole,
#   device = "png",
#   width = 8,
#   height = 6,
#   dpi = 300
# )