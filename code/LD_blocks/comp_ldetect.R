library(bigsnpr)
library(data.table)

ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr1.bed")
n_block0 <- nrow(ldetect_chr)


bim_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr1_filtered_200k.bim"
ukb_snps <- fread(bim_file, select = c(1, 2, 4), 
                  col.names = c("chr", "snp", "pos"))
setkey(ukb_snps, pos)
ldetect_chr = as.data.table(ldetect_chr)
ldetect_chr[, block_id := .I]

ukb_snps[, ldetect_block := NA_integer_]

for (i in 1:nrow(ldetect_chr)) {
  ukb_snps[
    pos >= ldetect_chr$start[i] & 
      pos < ldetect_chr$stop[i],
    ldetect_block := i
  ]
}

table(ukb_snps$ldetect_block, useNA = "always")

library(ggplot2)
library(data.table)
library(scales) # For formatting axes

# 1. Prepare LD block data
ldetect_chr[, block_size := stop - start]
ldetect_chr[, block_mid := (start + stop)/2]

# 2. Prepare SNP count data
snp_counts <- ukb_snps[, .(snp_count = sum(!is.na(ldetect_block))), by = ldetect_block]
ldetect_chr <- merge(ldetect_chr, snp_counts, by.x = "block_id", by.y = "ldetect_block", all.x = TRUE)

# 3. Create block size visualization
p_block_size <- ggplot(ldetect_chr, aes(x = block_mid, y = block_size)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
  geom_text(aes(label = block_id), vjust = -0.5, size = 3, color = "black") +
  scale_x_continuous(
    "Chromosome Position (bp)",
    labels = comma,
    breaks = ldetect_chr$block_mid,
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    "Block Size (bp)",
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(title = "LD Block Sizes on Chromosome 1") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# 4. Create SNP density visualization
p_snp_density <- ggplot(ldetect_chr, aes(x = block_mid, y = snp_count)) +
  geom_bar(stat = "identity", fill = "darkorange", width = 0.8) +
  geom_text(aes(label = block_id), vjust = -0.5, size = 3, color = "black") +
  scale_x_continuous(
    "Chromosome Position (bp)",
    labels = comma,
    breaks = ldetect_chr$block_mid,
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    "Number of SNPs",
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(title = "SNP Count per LD Block on Chromosome 1") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# 5. Create combined visualization
p_combined <- ggplot(ldetect_chr) +
  geom_bar(aes(x = block_mid, y = block_size), 
           stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_bar(aes(x = block_mid, y = snp_count * max(block_size)/max(snp_count)), 
           stat = "identity", fill = "darkorange", alpha = 0.7) +
  scale_x_continuous(
    "Chromosome Position (bp)",
    labels = comma,
    breaks = ldetect_chr$block_mid,
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    "Block Size (bp)",
    labels = comma,
    sec.axis = sec_axis(
      ~ . * max(snp_count)/max(block_size), 
      name = "Number of SNPs",
      labels = comma
    )
  ) +
  geom_text(aes(x = block_mid, y = block_size, label = block_id), 
            vjust = -0.5, size = 3, color = "black") +
  labs(title = "LD Block Sizes vs SNP Count on Chromosome 1") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y.right = element_text(color = "darkorange"),
    axis.text.y.right = element_text(color = "darkorange")
  )

# 6. Create dot plot showing size vs SNP count
p_correlation <- ggplot(ldetect_chr, aes(x = block_size, y = snp_count)) +
  geom_point(size = 3, color = "purple") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgray") +
  geom_text(aes(label = block_id), vjust = -1, size = 3) +
  scale_x_continuous("Block Size (bp)", labels = comma) +
  scale_y_continuous("Number of SNPs", labels = comma) +
  labs(title = "LD Block Size vs SNP Count Correlation") +
  theme_minimal()

# Print plots
print(p_block_size)
print(p_snp_density)
print(p_combined)
print(p_correlation)

# Save plots
ggsave(paste0(results_dir, "/ld_block_sizes.png"), 
       p_block_size, width = 12, height = 6)
ggsave(paste0(results_dir, "/snp_counts_per_block.png"), 
       p_snp_density, width = 12, height = 6)
ggsave(paste0(results_dir, "/combined_block_metrics.png"), 
       p_combined, width = 12, height = 6)
ggsave(paste0(results_dir, "/size_vs_snp_correlation.png"), 
       p_correlation, width = 8, height = 6)