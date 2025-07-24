library(data.table)
library(dplyr)
library(genio)
library(BEDMatrix)

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter" 
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

h2 <- 0.5
p <- 1000
n_snps <- 1000

# (1) Select random causal SNPs genome-wide --------------------------------
# Combine frequencies across chromosomes
all_freq <- rbindlist(lapply(1:1, function(chr) {
  fread(paste0(base_dir, "/ukb_chr", chr, "_filtered_200k.afreq"),
        select = c(1, 2, 3, 4, 6),  # CHR, SNP ID, REF, ALT, ALT_FREQ
        col.names = c("chr", "snp", "ref", "alt", "maf"))
}))
# Randomly select SNPs
set.seed(123)
causal_snps <- all_freq[sample(.N, n_snps)][, index := .I]
causal_snps_list <- causal_snps %>% select(snp)
# Write causal SNPs to extract
write.table(causal_snps_list, file = paste0(results_dir, "/causal_snps_chr1.txt"),
            col.names = FALSE, sep = "\t", quote = F, row.names = F)
# (2) Simulate effects -----------------------------------------------------
set.seed(123)
a <- rnorm(p, mean = 6, sd = sqrt(4))
d <- rep(0, length(a))

# (3) Process chromosomes ---------------------------------------------------
total_geno <- NULL
# Process chromosome 1
CHR <- 1
cat("Processing chromosome", CHR, "\n")
system(paste(
  plink2_path,
  "--bfile", paste0(base_dir, "/ukb_chr", CHR, "_filtered_200k"),
  "--extract", paste0(results_dir, "/causal_snps_chr1.txt"),
  "--make-bed",
  "--out", paste0(base_dir, "/temp_chr", CHR),
  "--memory 8000",
  "--threads 24"
))

plink_file <- paste0(base_dir, "/temp_chr", CHR)
bim_file <- paste0(plink_file, ".bim")
fam_file <- paste0(plink_file, ".fam")

# Read sample information
fam <- fread(fam_file, col.names = c("FID", "IID", "pat", "mat", "sex", "pheno"))

if(is.null(total_geno)) {
  total_geno <- data.frame(
    FID = fam$FID,
    IID = fam$IID,
    score = 0
  )
}

# Load genotype data with BEDMatrix (no transpose)
bed_path <- paste0(plink_file, ".bed")
geno <- BEDMatrix(bed_path, simple_names = TRUE)  # Individuals x SNPs

# Read BIM file
bim <- fread(bim_file, col.names = c("chr", "snp", "pos_cM", "pos", "alt", "ref"))

# Get causal SNPs for this chromosome
chr_snps <- causal_snps$snp  # All causal SNPs are on chr1

# Match genotype columns to causal SNPs
snps_in_bim <- bim$snp
geno_snps <- colnames(geno)
matched_snps <- intersect(chr_snps, geno_snps)

# Subset genotype matrix to matched SNPs
geno_sub <- as.matrix(geno[, matched_snps, drop = FALSE])

# Get effect sizes for matched SNPs
match_idx <- match(matched_snps, causal_snps$snp)
a_vals <- a[match_idx]
d_vals <- d[match_idx]

# Calculate genetic score using vectorized operations
for (j in 1:length(matched_snps)) {
  g <- geno_sub[, j]  # Genotypes for all individuals for this SNP
  g[is.na(g)] <- 1    # Set missing to heterozygous
  
  # Calculate contributions using vectorization
  contrib <- case_when(
    g == 0 ~ -a_vals[j],
    g == 1 ~ d_vals[j],
    g == 2 ~ a_vals[j],
    TRUE ~ 0
  )
  
  total_geno$score <- total_geno$score + contrib
}

# Clean up temporary files
unlink(paste0(base_dir, "/temp_chr", CHR, c(".bed", ".bim", ".fam", ".log")))

# (3) Final phenotype calculation --------------------------------------------
var_g <- var(total_geno$score)
var_e <- (1/h2 - 1) * var_g
total_geno$pheno <- total_geno$score + rnorm(nrow(total_geno), 0, sqrt(var_e))
pheno <- data.frame(
  FID = as.character(total_geno$FID),
  IID = as.character(total_geno$IID),
  pheno = total_geno$score + rnorm(nrow(total_geno), 0, sqrt((1/h2 - 1) * var(total_geno$score)))
)
# Save results with IDs
write.table(pheno,
            paste0(results_dir, "/simulated_pheno_chr1.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

system(paste(
  plink2_path,
  "--bfile", paste0(base_dir, "/ukb_chr", CHR, "_filtered_200k"),
  "--pheno", paste0(results_dir, "/simulated_pheno_chr1.txt"),
  "--glm no-x-sex allow-no-covars",
  "--out", paste0(results_dir, "/com_ldetect/chr", CHR, "_gwas_test"),
  "--memory 8000",
  "--threads 24"
))



###--------------PLOTs-------------------------
library(ggplot2)
# Read GWAS results
gwas_res = fread(paste0(results_dir, "/com_ldetect/chr1_gwas_test.pheno.glm.linear"))
# Prepare Manhattan plot data
manhattan <- gwas_res %>%
  rename(BP = POS, SNP = ID) %>%
  mutate(P = as.numeric(P)) %>%
  arrange(BP)  # Sort by position

# Identify causal SNPs (modify these based on your simulation)
causal_a_snps <- causal_snps$snp  # Additive causal SNPs
causal_d_snps <- c()              # Dominant causal SNPs (if any)

# Add causal SNP annotation
manhattan <- manhattan %>%
  mutate(is_causal = case_when(
    SNP %in% causal_d_snps ~ "dominant",
    SNP %in% causal_a_snps ~ "additive",
    TRUE ~ "non-causal"
  ))

# Read LD blocks
ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr1.bed")
# Assign colors to LD blocks
n_blocks <- nrow(ldetect_chr)
ldetect_chr$block_color <- ifelse(seq_len(n_blocks) %% 2 == 1, "gray90", "white")

# Create plot
p_gwas <- ggplot(manhattan, aes(x = BP, y = -log10(P))) +
  # Draw LD blocks as alternating rectangles
  geom_rect(
    data = ldetect_chr,
    aes(xmin = start, xmax = stop, ymin = -Inf, ymax = Inf, fill = block_color),
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +  # Use the colors directly from the data
  
  # Add all GWAS points
  geom_point(aes(color = is_causal), alpha = 0.8, size = 1.5) +
  
  # Color coding for causal SNPs
  scale_color_manual(
    name = "Causal SNPs",
    values = c("additive" = "orange", "dominant" = "green", "non-causal" = "black"),
    breaks = c("additive", "dominant")  # Only show these in legend
  ) +
  
  # Bonferroni correction line (genome-wide for this chromosome)
  geom_hline(
    yintercept = -log10(0.05/nrow(manhattan)),
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Axis labels and formatting
  labs(
    title = "Manhattan Plot - Chromosome 1",
    x = "Position (bp)",
    y = expression(-log[10](p-value))
  ) +
  scale_x_continuous(labels = scales::comma) +  # Format position labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Add space at top
  
  # Theme customization
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),  # 添加这行
    plot.background = element_rect(fill = "white", color = NA)   # 添加这行
  )

# Print plot
print(p_gwas)

# Save plot
ggsave(paste0(results_dir, "/manhattan_chr1.png"), 
       p_gwas, width = 12, height = 6, dpi = 300)
