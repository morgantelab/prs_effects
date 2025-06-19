library(data.table)
library(dplyr)
library(genio)

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter" 
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

h2 <- 0.5
p <- 1000
n_snps <- 1000

# (1) Select random causal SNPs genome-wide --------------------------------
# Combine frequencies across chromosomes
all_freq <- rbindlist(lapply(1:22, function(chr) {
  fread(paste0(base_dir, "/ukb_chr", chr, "_filtered_200k.afreq"),
        select = c(1, 2, 3, 4, 6),  # CHR, SNP ID, REF, ALT, ALT_FREQ
        col.names = c("chr", "snp", "ref", "alt", "maf"))
}))
# Randomly select SNPs
set.seed(123)
causal_snps <- all_freq[sample(.N, n_snps)][, index := .I]
causal_snps_list <- causal_snps %>% select(snp)
# Write causal SNPs to extract
write.table(causal_snps_list, file = paste0(results_dir, "/causal_snps_largesample_200k.txt"),
       col.names = FALSE, sep = "\t", quote = F, row.names = F)
# (2) Simulate effects -----------------------------------------------------
set.seed(123)
a <- rnorm(p, mean = 0, sd = sqrt(2))
set.seed(123)
dom_type <- rep(c("complete", "incomplete", "additive"), times = c(500, 300, 200))

d <- numeric(p)

for (i in 1:p) {
  if (dom_type[i] == "complete") {
    d[i] <- sample(c(1, -1), 1) * a[i]
  } else if (dom_type[i] == "incomplete") {
    d[i] <- runif(1, 0.2, 0.8) * sample(c(1, -1), 1) * a[i]
    } else if (dom_type[i] == "additive"){
    d[i] <- 0
  }
}

# (3) Process chromosomes ---------------------------------------------------
total_geno <- NULL
ID <- NULL

for(CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  system(paste(
    plink2_path,
    "--bfile", paste0(base_dir, "/ukb_chr", CHR, "_filtered_200k"),
    "--extract", paste0(results_dir, "/causal_snps_largesample_200k.txt"),
    "--make-bed",
    "--out", paste0(base_dir, "/temp_chr", CHR),
    "--memory 8000",
    "--threads 8"
  ))

  plink_data <- read_plink(paste0(base_dir, "/temp_chr", CHR))

  # Get sample IDs (FID and IID)
  if(is.null(ID)) {
    ID <- plink_data$fam[, c("fam", "id")]
    colnames(ID) <- c("FID", "IID")
    total_geno <- data.frame(
      FID = ID$FID,
      IID = ID$IID,
      score = 0
    )
  }

  geno_matrix <- t(plink_data$X)
  colnames(geno_matrix) <- plink_data$bim$id

  chr_snps <- causal_snps[chr == CHR, snp]
  geno_matrix <- geno_matrix[, chr_snps, drop = FALSE]

  a_vals <- a[causal_snps[chr == CHR, index]]
  d_vals <- d[causal_snps[chr == CHR, index]]

  for (j in 1:ncol(geno_matrix)) {
    for (i in 1:nrow(geno_matrix)) {
      genotype <- geno_matrix[i, j]
      if (is.na(genotype)) next

      if (geno_matrix[i, j] == 0) {
        total_geno$score[i] <- total_geno$score[i] - a_vals[j]
      } else if (geno_matrix[i, j] == 1) {
        total_geno$score[i] <- total_geno$score[i] + d_vals[j]
      } else if (geno_matrix[i, j] == 2) {
        total_geno$score[i] <- total_geno$score[i] + a_vals[j]
      }
    }
  }

  unlink(paste0(base_dir, "/temp_chr", CHR, c(".bed", ".bim", ".fam", ".log")))
}

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
            paste0(results_dir, "/simulated_pheno_large_200k.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(total_geno[, c("FID", "IID", "score")],
            paste0(results_dir, "/simulated_geno_large_200k.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("/n Success! Simulated pheno and geno saved to:", results_dir, "\n")



  
### 3. Calculate PRS ------

# calculate_prs <- function(geno_data, effect_results, alpha=0.05, total_loci) {
#   # 1. Identify SNPs significant for EITHER effect (Bonferroni corrected)
#   sig_idx <- which(
#     (effect_results$p_a < (alpha/total_loci)) |
#       (effect_results$p_d < (alpha/total_loci))
#   )
#   sig_SNPs <- effect_results$SNP[sig_idx]
# 
#   # 2. Initialize PRS vector
#   prs <- numeric(nrow(geno_data))
# 
#   # 3. Calculate PRS using explicit genotype conditions
#   for(i in seq_along(sig_SNPs)) {
#     snp <- sig_SNPs[i]
#     genotypes <- geno_data[[snp]]
#     a <- effect_results$a_hat[sig_idx[i]]
#     d <- effect_results$d_hat[sig_idx[i]]
# 
#     for(j in 1:length(genotypes)) {
#       gt <- genotypes[j]
# 
#       if (gt == 2) {                # Homozygous alternate (AA)
#         prs[j] <- prs[j] + a
#       } else if (gt == 1) {         # Heterozygous (Aa)
#         prs[j] <- prs[j] + d
#       } else if (gt == 0) {         # Homozygous reference (aa)
#         prs[j] <- prs[j] - a
#       }
#     }
#   }
# 
#   return(prs)
# }
# 
# 
# prs_effects <- calculate_prs(
#   geno_data = full_data[, -c("FID", "IID", "pheno")],  # Remove non-genotype columns
#   effect_results = effect_results,
#   alpha = 0.05,
#   total_loci = nrow(effect_results)
# )
# 
# 
# final_results <- data.table(
#   FID = full_data$FID,
#   IID = full_data$IID,
#   Phenotype = full_data$pheno,
#   PRS_Effect = prs_effects
# )
# 
# output_file <- paste0(results_dir, "/prs_effects.csv")
# write.table(final_results,
#             file = output_file,
#             row.names = FALSE,
#             col.names = TRUE,
#             quote = FALSE,
#             sep = "\t")
# 
# cat("Analysis completed successfully! Results saved to prs_results.csv\n")


