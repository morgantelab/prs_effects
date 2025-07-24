library(data.table)
library(ggplot2)
library(bigsnpr)

# Set paths
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/LD_decay"
geno_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/filter"
# Load LD blocks
ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr1.bed")

output_prefix_50kb <-paste0(results_dir, "/chr1_10kb_200k_maf0.1")
bfile_path <- paste0(geno_dir, "/ukb_chr1_filtered")
keep_file_path<-"/data2/morgante_lab/ukbiobank_projects/prs_effects/results/ID_filtered.txt"
plink_cmd <- sprintf(
  "%s --bfile %s --chr 1 --r2-unphased --maf 0.1 --geno 0.1 --hwe 10E-5 --mac 5 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 10 --keep %s --out %s",
  plink2_path,
  bfile_path,
  keep_file_path,
  output_prefix_50kb
)

system(plink_cmd)

library(tidyverse)
file_in <- paste0(output_prefix_50kb, ".vcor")
data <- fread(file_in, header = T) %>%
  dplyr::select(CHR = "#CHROM_A", POS_A, POS_B, UNPHASED_R2)

colnames(data) <- c("chr", "snp_1", "snp_2", "r2")
data$distance <- abs(data$snp_1 - data$snp_2)
data$distance_kb <- abs(data$snp_1 - data$snp_2) / 1000
# Create bins
data$bin_dist <- ggplot2::cut_interval(data$distance, n = 500)

# Get left boundary
data$bin <- ggplot2::cut_interval(data$distance, n = 500)
data$bin_bdr <- as.numeric(stringr::str_split(data$bin, ",", simplify = T)[,1 ] %>%
                             stringr::str_replace_all("\\(", "") %>%
                             stringr::str_replace_all("\\[", ""))

# Group by bin and get mean R^2
out <- data %>%
  dplyr::group_by(bin_bdr) %>%
  dplyr::summarise(mean = mean(r2, na.rm = T))

# Write table
file_out <- file.path(results_dir, paste0("/chr1_10kb.txt"))
write.table(out,
            file = file_out,
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")

out$bin_bdr_kb <- out$bin_bdr / 1000  

r2_plot_main = out %>% ggplot() +
  geom_line(aes(bin_bdr_kb, mean)) +
  xlab("Distance between SNPs (kb)") +
  ylab(bquote(.("Mean r")^2)) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, .8)) +
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     limits = c(0, 0.6))
ggsave(plot = r2_plot_main, filename = paste0(results_dir, "/chr1_10kb_200k_maf0.1.png"))