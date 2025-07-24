library(dplyr)
library(ggplot2)
library(data.table)
library(readxl)

# (6) Plots ------------------------------------------------------------

### 1. Manhattan Plot
effect_results = fread(paste0(results_dir,"/effect_size_large_200k.txt"))

causal_d_snps <- causal_snps$snp[dom_type[causal_snps$index] %in% c("complete", "incomplete")]
causal_a_snps <- causal_snps$snp

manhattan_data <- effect_results %>%
  filter(chromosome == 1) %>%
  mutate(is_causal = case_when(
    SNP %in% causal_d_snps ~ "dominant",
    SNP %in% causal_a_snps & !SNP %in% causal_d_snps ~ "additive",
    TRUE ~ "non-causal"  # All other SNPs are labeled as non-causal
  ))
manhattan_data_a <- manhattan_data %>%
  rename(P = p_a) %>%
  rename(CHR = chr) %>%
  rename(BP = position)
# png(file.path(results_dir, "manhattan_additive.png"), width=12, height=6, units="in", res=300)
#
# manhattan(
#   manhattan_data_a,
#   main = "Additive Effects",
#   suggestiveline = -log10(1e-5),
#   genomewideline = -log10(0.05/nrow(manhattan_data)),
#   highlight = causal_a_snps,
#   col = c("black", "grey")
# )
# dev.off()

manhattan_data_d <- manhattan_data %>%
  rename(P = p_d) %>%
  rename(CHR = chr) %>%
  rename(BP = position)
# png(file.path(results_dir, "manhattan_dominant.png"), width=12, height=6, units="in", res=300)

# manhattan(
#   manhattan_data_d,
#   main = "Domiance Effects",
#   suggestiveline = -log10(1e-5),
#   genomewideline = -log10(0.05/nrow(manhattan_data)),
#   highlight = causal_d_snps,
#   col = c("black", "grey")
# )
# dev.off()

library(ggplot2)
additive_data <- manhattan_data_a %>%

  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(manhattan_data_a, ., by=c("CHR"="CHR")) %>%

  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = additive_data %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p_add <- ggplot(additive_data, aes(x=BPcum, y=-log10(P))) +

  # Show all points
  geom_point(color = "gray70", alpha = 0.8, size = 2) +

  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

  # Add highlighted points
  geom_point(data = subset(additive_data, !is.na(is_causal)),aes(color = is_causal),size = 2) +
  scale_color_manual( name = "Causal SNPs",
                      values = c("additive" = "green", "dominant" = "orange")) +
  # Add Bonferroni correction significant line
  geom_hline(yintercept = -log10(0.05/nrow(all_freq)),
             color = "red", linetype = "dashed") +

  # Add labs
  labs(
    title = "Additive Effects",
    x = "Chromosome",
    y = expression(-log[10](p-value))) +
  # Custom the theme:
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )


dominance_data <- manhattan_data_d %>%

  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(manhattan_data_d, ., by=c("CHR"="CHR")) %>%

  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = dominance_data %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p_dom <- ggplot(dominance_data, aes(x=BPcum, y=-log10(P))) +

  # Show all points
  geom_point(color = "gray70", alpha = 0.8, size = 2) +

  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

  # Add highlighted points
  geom_point(data = subset(dominance_data, !is.na(is_causal)),aes(color = is_causal),size = 2) +
  scale_color_manual( name = "Causal SNPs",
                      values = c("additive" = "green", "dominant" = "orange")) +
  # Add Bonferroni correction significant line
  geom_hline(yintercept = -log10(0.05/nrow(dominance_data)),
             color = "red", linetype = "dashed") +

  # Add labs
  labs(
    title = "Dominance Effects",
    x = "Chromosome",
    y = expression(-log[10](p-value))) +
  # Custom the theme:
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )

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

ld_blocks <- ldetect_chr %>%
  mutate(
    start_cum = start + additive_data$tot[1],  # Cumulative adjustment (0 for chr1)
    end_cum = stop + additive_data$tot[1]
  )

p_add <- ggplot(additive_data, aes(x = BPcum, y = -log10(P))) +
  # Add LD blocks first (background layer)
  geom_rect(
    data = ld_blocks,
    aes(xmin = start_cum, xmax = end_cum, ymin = -Inf, ymax = Inf),
    fill = ld_blocks$block_color,
    inherit.aes = FALSE
  ) +
  # Original plot elements
  geom_point(color = "gray70", alpha = 0.8, size = 2) +
  geom_point(data = subset(additive_data, !is.na(is_causal)), 
             aes(color = is_causal), size = 2) +
  scale_color_manual( name = "Causal SNPs",
                      values = c("additive" = "green", "dominant" = "orange")) +
  # Add Bonferroni correction significant line
  geom_hline(yintercept = -log10(0.05/nrow(all_freq)),
             color = "red", linetype = "dashed") +
  
  # Add labs
  labs(
    title = "Additive Effects",
    x = "Chromosome",
    y = expression(-log[10](p-value))) +
  # Custom the theme:
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )
