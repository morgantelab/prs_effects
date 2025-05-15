library(data.table)
library(dplyr)
library(genio)

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter" 
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"

# h2 <- 0.5
# p <- 1000
# n_snps <- 1000
# 
# set.seed(123)  # For reproducibility
# full_fam <- fread(
#  paste0(base_dir, "/ukb_chr21_filtered.fam"),
#   col.names = c("FID", "IID", "PID", "MID", "Sex", "Pheno")
# )
# 
# sample_indices <- sample(nrow(full_fam), 50000)
# sample_ids <- full_fam[sample_indices, .(FID, IID)]
# # Write sample IDs to keep
# fwrite(sample_ids[, .(FID, IID)],
#        paste0(results_dir, "/subset_samples.txt"),
#        sep = " ", col.names = FALSE)
# 
# # (1) Select random causal SNPs genome-wide --------------------------------
# # Combine frequencies across chromosomes
# all_freq <- rbindlist(lapply(1:22, function(chr) {
#   fread(paste0(base_dir, "/ukb_chr", chr, "_filtered.afreq"),
#         select = c(1, 2, 3, 4, 6),  # CHR, SNP ID, REF, ALT, ALT_FREQ
#         col.names = c("chr", "snp", "ref", "alt", "maf"))
# }))
# ### 339047 SNPs left
# # Randomly select SNPs
# set.seed(123)
# causal_snps <- all_freq[sample(.N, n_snps)][, index := .I]
# causal_snps_list <- causal_snps %>% select(snp)
# # Write causal SNPs to extract
# write.table(causal_snps_list, file = paste0(results_dir, "/causal_snps.txt"),
#        col.names = FALSE, sep = "\t", quote = F, row.names = F)
# # (2) Simulate effects -----------------------------------------------------
# set.seed(123)
# a <- rnorm(p, mean = 0, sd = sqrt(2))
# set.seed(123)
# dom_type <- rep(c("complete", "incomplete", "additive"), times = c(500, 300, 200))
# 
# d <- numeric(p)
# 
# for (i in 1:p) {
#   if (dom_type[i] == "complete") {
#     d[i] <- sample(c(1, -1), 1) * a[i]
#   } else if (dom_type[i] == "incomplete") {
#     d[i] <- runif(1, 0.2, 0.8) * sample(c(1, -1), 1) * a[i]
#     } else if (dom_type[i] == "additive"){
#     d[i] <- 0
#   }
# }
# 
# # (3) Process chromosomes ---------------------------------------------------
# total_geno <- NULL
# ID <- NULL
# 
# for(CHR in 1:22) {
#   cat("Processing chromosome", CHR, "\n")
#   system(paste(
#     plink2_path,
#     "--bfile", paste0(base_dir, "/ukb_chr", CHR, "_filtered"),
#     "--keep", paste0(results_dir, "/subset_samples.txt"),
#     "--extract", paste0(results_dir, "/causal_snps.txt"),
#     "--make-bed",
#     "--out", paste0(base_dir, "/temp_chr", CHR),
#     "--memory 8000",
#     "--threads 8"
#   ))
# 
#   plink_data <- read_plink(paste0(base_dir, "/temp_chr", CHR))
# 
#   # Get sample IDs (FID and IID)
#   if(is.null(ID)) {
#     ID <- plink_data$fam[, c("fam", "id")]  # Use correct column names from fam file
#     colnames(ID) <- c("FID", "IID")
#     total_geno <- data.frame(
#       FID = ID$FID,
#       IID = ID$IID,
#       score = 0
#     )
#   }
# 
#   geno_matrix <- t(plink_data$X)
#   colnames(geno_matrix) <- plink_data$bim$id
# 
#   chr_snps <- causal_snps[chr == CHR, snp]
#   geno_matrix <- geno_matrix[, chr_snps, drop = FALSE]
# 
#   a_vals <- a[causal_snps[chr == CHR, index]]
#   d_vals <- d[causal_snps[chr == CHR, index]]
# 
#   for (j in 1:ncol(geno_matrix)) {
#     for (i in 1:nrow(geno_matrix)) {
#       genotype <- geno_matrix[i, j]
#       if (is.na(genotype)) next
# 
#       if (geno_matrix[i, j] == 0) {
#         total_geno$score[i] <- total_geno$score[i] - a_vals[j]
#       } else if (geno_matrix[i, j] == 1) {
#         total_geno$score[i] <- total_geno$score[i] + d_vals[j]
#       } else if (geno_matrix[i, j] == 2) {
#         total_geno$score[i] <- total_geno$score[i] + a_vals[j]
#       }
#     }
#   }
# 
#   unlink(paste0(base_dir, "/temp_chr", CHR, c(".bed", ".bim", ".fam", ".log")))
# }
# 
# # (3) Final phenotype calculation --------------------------------------------
# var_g <- var(total_geno$score)
# var_e <- (1/h2 - 1) * var_g
# total_geno$pheno <- total_geno$score + rnorm(nrow(total_geno), 0, sqrt(var_e))
# pheno <- data.frame(
#   FID = as.character(total_geno$FID),
#   IID = as.character(total_geno$IID),
#   pheno = total_geno$score + rnorm(nrow(total_geno), 0, sqrt((1/h2 - 1) * var(total_geno$score)))
# )
# # Save results with IDs
# write.table(pheno,
#             paste0(results_dir, "/simulated_pheno_small.txt"),
#             row.names = FALSE, quote = FALSE, sep = "\t")
# 
# write.table(total_geno[, c("FID", "IID", "score")],
#             paste0(results_dir, "/simulated_geno_small.txt"),
#             row.names = FALSE, quote = FALSE, sep = "\t")
# 
# cat("Success! Results saved to:", results_dir, "\n")


###-----------------------Directly Read in Pheno, All Above can not run-----------------------------------------
pheno <- fread(paste0(results_dir, "/simulated_pheno_small.txt"), 
               colClasses = c("character", "character", "numeric"))  

# (4) Estimating Effects -----------------------------------------------------

### 1. Used Functions

calculate_effects <- function(data, locus_col) {
  # Convert genotypes to ordered factor
  data[[locus_col]] <- factor(data[[locus_col]], 
                              levels = c("0", "1", "2"),
                              ordered = TRUE)
  
  # Calculate means with NA removal
  mean_genotypes <- tapply(data$pheno, data[[locus_col]], mean, na.rm = TRUE)
  
  # Calculate effects
  a <- 0.5 * (mean_genotypes["2"] - mean_genotypes["0"])
  d <- mean_genotypes["1"] - 0.5 * (mean_genotypes["2"] + mean_genotypes["0"])
  
  return(data.frame(SNP = locus_col, a = a, d = d))
}


calculate_welch <- function(data, snp_cols) {
  n_snps <- length(snp_cols)
  results <- data.frame(
    SNP = snp_cols,
    a_hat = numeric(n_snps),
    d_hat = numeric(n_snps),
    p_a = numeric(n_snps),
    p_d = numeric(n_snps)
  )
  
  for(j in 1:n_snps) {
    snp <- snp_cols[j]
    X <- data[[snp]]
    pheno <- data$pheno
    
    AA <- X == 2
    Aa <- X == 1
    aa <- X == 0
    
    n_AA <- sum(AA, na.rm = TRUE)
    n_Aa <- sum(Aa, na.rm = TRUE)
    n_aa <- sum(aa, na.rm = TRUE)
    
    # Additive effect calculations
    y_AA <- pheno[AA]
    y_aa <- pheno[aa]
    
    mean_AA <- mean(y_AA, na.rm = TRUE)
    mean_aa <- mean(y_aa, na.rm = TRUE)
    var_AA <- var(y_AA, na.rm = TRUE)
    var_aa <- var(y_aa, na.rm = TRUE)
    
    a_hat <- (mean_AA - mean_aa)/2
    var_a <- 0.25*(var_AA/n_AA + var_aa/n_aa)
    df_a <- (var_AA/n_AA + var_aa/n_aa)^2 / 
      ((var_AA/n_AA)^2/(n_AA-1) + (var_aa/n_aa)^2/(n_aa-1))
    t_a <- a_hat/sqrt(var_a)
    p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)
    
    # Dominance effect calculations
    y_Aa <- pheno[Aa]
    mean_Aa <- mean(y_Aa, na.rm = TRUE)
    var_Aa <- var(y_Aa, na.rm = TRUE)
    
    d_hat <- mean_Aa - (mean_AA + mean_aa)/2
    var_d <- var_Aa/n_Aa + 0.25*(var_AA/n_AA + var_aa/n_aa)
    df_d <- (var_d)^2 / (
      (var_Aa/n_Aa)^2/(n_Aa-1) + 
        (var_AA/(4*n_AA))^2/(n_AA-1) + 
        (var_aa/(4*n_aa))^2/(n_aa-1)
    )
    t_d <- d_hat/sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}



### 2. Loop over chromosome
process_geno_files <- function(base_dir, pheno) {
  # Initialize list to store results
  all_results <- list()
  
  for(chr in 1:22) {
    cat("\nProcessing chromosome", chr, "\n")
    
    # Read PLINK files
    plink_data <- read_plink(paste0(base_dir, "/ukb_chr", chr, "_filtered"))
    cat("Loaded", nrow(plink_data$fam), "samples and", nrow(plink_data$bim), "SNPs \n")
    
    # Create genotype data frame
    geno_df <- as.data.frame(t(plink_data$X))
    colnames(geno_df) <- plink_data$bim$id  # SNP IDs as column names
    cat("\nTotal of", ncol(geno_df), "SNPs are included", "\n")
    geno_df <- cbind(plink_data$fam[c("fam", "id")], geno_df)
    names(geno_df)[1:2] <- c("FID", "IID")
    
    snp_info <- plink_data$bim %>%
      dplyr::select(id, chr, pos, allele1, allele2) %>%
      dplyr::rename(SNP = id, 
                    chromosome = chr,
                    position = pos,
                    ref_allele = allele1,
                    alt_allele = allele2)
    
    # Merge with phenotype data
    full_data <- merge(pheno, geno_df, by = c("FID", "IID"))
    
    # Get SNP columns (excluding ID and phenotype columns)
    snp_cols <- setdiff(colnames(full_data), c("FID", "IID", "pheno"))
    
    # Calculate Welch statistics
    chr_results <- calculate_welch(full_data, snp_cols)
    chr_results$chr <- chr
    
    chr_results <- chr_results %>%
      left_join(snp_info, by = "SNP") %>%
      dplyr::select(chromosome, position, SNP, ref_allele, alt_allele, everything())
    
    all_results[[as.character(chr)]] <- chr_results
    
    rm(plink_data, geno_df, full_data)
    gc()
    
      }
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  return(combined_results)
}
start_time <- Sys.time()

cat("\nStarting estimation at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
effect_results = process_geno_files(base_dir, pheno)
end_time <- Sys.time()

cat("\nCompleted estimation at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

write.table(effect_results, paste0(results_dir, "/effect_size.txt"),
                           row.names = FALSE, quote = FALSE, sep = "\t")

cat("Success! Results saved to:", results_dir, "\n")
  
# ### 3. Calculate PRS ------
# 
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
# 
# 

# # (6) Plots and Matrix ------------------------------------------------------------
# 
# ### 1. Manhattan Plot
# 
# effect_results = fread(paste0(results_dir,"/effect_size.txt"))
# effect_results = merge(effect_results, all_freq[c("chr", )])
# library("qqman")
# 
# # Identify causal SNPs
# d_causal_indices <- which(dom_type %in% c("complete", "incomplete"))  # First 800 SNPs (500+300)
# 
# manhattan_a <- effect_results[, .(SNP, CHR = chr, BP = pos, P = p_a,
#                                   is_causal = SNP %in% causal_snps$snp)]
# 
# manhattan_d <- effect_results[, .(SNP, CHR = chr, BP = pos, P = p_d,
#                                   is_causal = .I %in% d_causal_indices)]
# 
# output_file <- file.path(results_dir, "manhattan_plots.pdf")
# 
# pdf(output_file, width = 12, height = 6)
# 
# par(mfrow = c(1, 2))
# 
# # Additive effects (all p causal SNPs)
# manhattan(manhattan_a, main = "Additive Effects",
#           genomewideline = -log10(0.05/nrow(manhattan_a)),
#           col = c("#276FBF", "#183059"))
# points(manhattan_a[is_causal == TRUE, .(BP, P)],
#        col = "red", pch = 1, cex = 1.2, lwd = 1.5)
# 
# # Dominance effects (only complete/incomplete dominance SNPs)
# manhattan(manhattan_d, main = "Dominance Effects",
#           genomewideline = -log10(0.05/nrow(manhattan_d)),
#           col = c("#276FBF", "#183059"))
# points(manhattan_d[is_causal == TRUE, .(BP, P)],
#        col = "red", pch = 1, cex = 1.2, lwd = 1.5)
# 
# # Add legend
# legend("topright",
#        legend = c("Causal SNPs", "Bonferroni Threshold"),
#        col = c("red", "black"),
#        pch = c(1, NA),
#        lty = c(NA, 2))
# dev.off()
# 
# message("Manhattan plots saved as 'manhattan_plots.pdf' in your results directory")
# 
# output_file <- paste0(results_dir, "/effects_dat.csv")
# write.table(effect_results,
#             file = output_file,
#             row.names = FALSE,
#             col.names = TRUE,
#             quote = FALSE,
#             sep = "\t")

