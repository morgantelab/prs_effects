library(data.table)
library(dplyr)
library(genio)

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter" 
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
h2 <- 0.5
p <- 1000
n_snps <- 1000

###-----------------------Directly Read in Pheno, All Above do not need to run-----------------------------------------
pheno <- fread(paste0(results_dir, "/simulated_pheno_large_200k.txt"),
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
    plink_data <- read_plink(paste0(base_dir, "/ukb_chr", chr, "_filtered_200k"))
    cat("Loaded", nrow(plink_data$fam), "samples and", nrow(plink_data$bim), "SNPs \n")
    
    # Create genotype data frame
    geno_df <- as.data.frame(t(plink_data$X))
    colnames(geno_df) <- plink_data$bim$id  # SNP IDs as column names
    cat("\nTotal of", ncol(geno_df), "SNPs are included", "\n")
    geno_df <- cbind(plink_data$fam[c("fam", "id")], geno_df)
    names(geno_df)[1:2] <- c("FID", "IID")
    
    snp_info <- plink_data$bim %>%
      dplyr::select(id, chr, pos, ref, alt) %>%
      dplyr::rename(SNP = id, 
                    chromosome = chr,
                    position = pos,
                    ref_allele = ref,
                    alt_allele = alt)
    
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

write.table(effect_results, paste0(results_dir, "/effect_size_large_200k.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("Success! Effecet Size results saved to:", results_dir, "\n")


# # (6) Plots ------------------------------------------------------------
# 
# ### 1. Manhattan Plot
# 
# effect_results = fread(paste0(results_dir,"/effect_size_large_200k.txt"))
# 
# causal_d_snps <- causal_snps$snp[dom_type[causal_snps$index] %in% c("complete", "incomplete")]
# causal_a_snps <- causal_snps$snp
# 
# # library("qqman")
# 
# # Prepare data for qqman
# manhattan_data <- effect_results %>%
#   mutate(is_causal = case_when(
#     SNP %in% causal_d_snps ~ "dominant",
#     SNP %in% causal_a_snps & !SNP %in% causal_d_snps ~ "additive",
#     TRUE ~ "non-causal"  # All other SNPs are labeled as non-causal
#   ))
# manhattan_data_a <- manhattan_data %>% 
#   rename(P = p_a) %>%
#   rename(CHR = chr) %>%
#   rename(BP = position)
# # png(file.path(results_dir, "manhattan_additive.png"), width=12, height=6, units="in", res=300)
# # 
# # manhattan(
# #   manhattan_data_a,
# #   main = "Additive Effects",
# #   suggestiveline = -log10(1e-5),
# #   genomewideline = -log10(0.05/nrow(manhattan_data)),
# #   highlight = causal_a_snps,
# #   col = c("black", "grey")
# # )
# # dev.off()
# 
# manhattan_data_d <- manhattan_data %>% 
#   rename(P = p_d) %>%
#   rename(CHR = chr) %>%
#   rename(BP = position)
# # png(file.path(results_dir, "manhattan_dominant.png"), width=12, height=6, units="in", res=300)
# 
# # manhattan(
# #   manhattan_data_d,
# #   main = "Domiance Effects",
# #   suggestiveline = -log10(1e-5),
# #   genomewideline = -log10(0.05/nrow(manhattan_data)),
# #   highlight = causal_d_snps,
# #   col = c("black", "grey")
# # )
# # dev.off()
# 
# library(ggplot2)
# additive_data <- manhattan_data_a %>% 
#   
#   # Compute chromosome size
#   group_by(CHR) %>% 
#   summarise(chr_len=max(BP)) %>% 
#   
#   # Calculate cumulative position of each chromosome
#   mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
#   select(-chr_len) %>%
#   
#   # Add this info to the initial dataset
#   left_join(manhattan_data_a, ., by=c("CHR"="CHR")) %>%
#   
#   # Add a cumulative position of each SNP
#   arrange(CHR, BP) %>%
#   mutate( BPcum=BP+tot)
# 
# axisdf = additive_data %>%
#   group_by(CHR) %>%
#   summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# p_add <- ggplot(additive_data, aes(x=BPcum, y=-log10(P))) +
#   
#   # Show all points
#   geom_point(color = "gray70", alpha = 0.8, size = 2) +
#   
#   # custom X axis:
#   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
#   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
#   
#   # Add highlighted points
#   geom_point(data = subset(additive_data, !is.na(is_causal)),aes(color = is_causal),size = 2) +
#   scale_color_manual( name = "Causal SNPs",
#                       values = c("additive" = "green", "dominant" = "orange")) +
#   # Add Bonferroni correction significant line
#   geom_hline(yintercept = -log10(0.05/nrow(all_freq)), 
#              color = "red", linetype = "dashed") +
#   
#   # Add labs
#   labs(
#     title = "Additive Effects",
#     x = "Chromosome",
#     y = expression(-log[10](p-value))) +
#   # Custom the theme:
#   theme(
#     legend.position = "right", 
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# ggsave(file.path(results_dir, "manhattan_additive_large.png"), p_add, width=12, height=6)
# 
# 
# dominance_data <- manhattan_data_d %>% 
#   
#   # Compute chromosome size
#   group_by(CHR) %>% 
#   summarise(chr_len=max(BP)) %>% 
#   
#   # Calculate cumulative position of each chromosome
#   mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
#   select(-chr_len) %>%
#   
#   # Add this info to the initial dataset
#   left_join(manhattan_data_d, ., by=c("CHR"="CHR")) %>%
#   
#   # Add a cumulative position of each SNP
#   arrange(CHR, BP) %>%
#   mutate( BPcum=BP+tot)
# 
# axisdf = dominance_data %>%
#   group_by(CHR) %>%
#   summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# p_dom <- ggplot(dominance_data, aes(x=BPcum, y=-log10(P))) +
#   
#   # Show all points
#   geom_point(color = "gray70", alpha = 0.8, size = 2) +
#   
#   # custom X axis:
#   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
#   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
#   
#   # Add highlighted points
#   geom_point(data = subset(dominance_data, !is.na(is_causal)),aes(color = is_causal),size = 2) +
#   scale_color_manual( name = "Causal SNPs",
#                       values = c("additive" = "green", "dominant" = "orange")) +
#   # Add Bonferroni correction significant line
#   geom_hline(yintercept = -log10(0.05/nrow(dominance_data)), 
#              color = "red", linetype = "dashed") +
#   
#   # Add labs
#   labs(
#     title = "Dominance Effects",
#     x = "Chromosome",
#     y = expression(-log[10](p-value))) +
#   # Custom the theme:
#   theme(
#     legend.position = "right", 
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(), 
#     axis.line = element_line(colour = "black"),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# ggsave(file.path(results_dir, "manhattan_domiance_large.png"), p_dom, width=12, height=6)
# 
# 
# # (7) Matrix ------------------------------------------------------------
# alpha <- 0.05/nrow(all_freq)
# effect_results <- effect_results %>% 
#   mutate(causal_a = ifelse(SNP %in% causal_a_snps, TRUE, FALSE)) %>%
#   mutate(causal_d = ifelse(SNP %in% causal_d_snps, TRUE, FALSE))
# 
# TP_a <- sum(effect_results$p_a[effect_results$causal_a] < alpha)       # True causal called significant
# FP_a <- sum(effect_results$p_a[!effect_results$causal_a] < alpha)      # Non-causal called significant
# FN_a <- sum(effect_results$p_a[effect_results$causal_a] >= alpha)      # True causal NOT called significant
# TN_a <- sum(effect_results$p_a[!effect_results$causal_a] >= alpha)     # Non-causal NOT called significant
# 
# 
# 
# TP_d <- sum(effect_results$p_d[effect_results$causal_d] < alpha)       # True causal called significant
# FP_d <- sum(effect_results$p_d[!effect_results$causal_d] < alpha)      # Non-causal called significant
# FN_d <- sum(effect_results$p_d[effect_results$causal_d] >= alpha)      # True causal NOT called significant
# TN_d <- sum(effect_results$p_d[!effect_results$causal_d] >= alpha)     # Non-causal NOT called significant
# 
# 
# 
# calculate_metrics <- function(TP, FP, FN, TN) {
#   data.frame(
#     TPR = TP/(TP + FN),
#     FDR = FP/(TP + FP),
#     TNR = TN/(TN + FP),
#     Precision = TP/(TP + FP),
#     Accuracy = (TP + TN)/(TP + FP + FN + TN),
#     Power = TP/(TP + FN)
#   )
# }
# 
# # Create matrices
# additive_matrix <- matrix(c(TP_a, FP_a, FN_a, TN_a), ncol = 2,
#                           dimnames = list(c("Causal", "Non-causal"),
#                                           c("significant", "Not significant")))
# 
# dominance_matrix <- matrix(c(TP_d, FP_d, FN_d, TN_d), ncol = 2,
#                            dimnames = list(c("Causal", "Non-causal"),
#                                            c("significant", "Not significant")))
# 
# # Calculate metrics
# additive_metrics <- calculate_metrics(TP_a, FP_a, FN_a, TN_a)
# dominance_metrics <- calculate_metrics(TP_d, FP_d, FN_d, TN_d)
# 
# cat("\nCombined Results:\n")
# list(
#   Additive = list(
#     Confusion_Matrix = matrix(c(TP_a, FP_a, FN_a, TN_a), ncol = 2,
#                               dimnames = list(c("Causal", "Non-causal"),
#                                               c("Significant", "Not significant"))),
#     Metrics = additive_metrics
#   ),
#   Dominance = list(
#     Confusion_Matrix = matrix(c(TP_d, FP_d, FN_d, TN_d), ncol = 2,
#                               dimnames = list(c("Causal", "Non-causal"),
#                                               c("Significant", "Not significant"))),
#     Metrics = dominance_metrics
#   )
# )
