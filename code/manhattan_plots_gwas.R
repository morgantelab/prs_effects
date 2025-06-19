library(data.table)
library(dplyr)
library(ggplot2)

results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/gwas_results/"

gwas_list <- vector("list", 22)

for (chr in 1:22) {
  chr_file <- file.path(results_dir, sprintf("chr%d_gwas_large.pheno.glm.linear", chr))
  
  if (file.exists(chr_file)) {
    cat(sprintf("Reading chromosome %d...\n", chr))
    gwas_list[[chr]] <- fread(chr_file)
    
    # Add chromosome number if not present
    if (!"CHR" %in% names(gwas_list[[chr]])) {
      gwas_list[[chr]][, CHR := chr]
    }
    
  }
}

# Combine all results (rbindlist is faster than do.call(rbind))
combined_results <- rbindlist(gwas_list, use.names = TRUE, fill = TRUE)

manhattan<- combined_results %>%
  rename(BP = POS) %>%
  rename(SNP = ID) %>%
  mutate(P = as.numeric(P))

manhattan <- manhattan%>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len = max(BP)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(manhattan, ., by = "CHR") %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  mutate(is_causal = case_when(
    SNP %in% causal_d_snps ~ "dominant",
    SNP %in% causal_a_snps & !SNP %in% causal_d_snps ~ "additive",
    TRUE ~ "non-causal"  # All other SNPs are labeled as non-causal
  ))


axisdf = manhattan %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p_gwas <- ggplot(manhattan, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point(color = "gray70", alpha = 0.8, size = 2) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data = subset(manhattan, !is.na(is_causal)),aes(color = is_causal),size = 2) +
  scale_color_manual(
    name = "Causal SNPs",
    values = c("additive" = "green", "dominant" = "orange")
  ) + 
  # Add Bonferroni correction significant line
  geom_hline(yintercept = -log10(0.05/nrow(manhattan)), 
             color = "red", linetype = "dashed") +
  
  # Add labs
  labs(
    title = "GWAS Results",
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

ggsave(file.path(results_dir, "manhattan_gwas_large.png"), p_gwas, width=12, height=6)


# alpha <- 0.05/nrow(combined_results)
# combined_results <- combined_results %>% 
#   mutate(causal_a = ifelse(ID %in% causal_a_snps, TRUE, FALSE)) %>%
#   mutate(causal_d = ifelse(ID %in% causal_d_snps, TRUE, FALSE))
# 
# TP_a <- sum(combined_results$p[combined_results$causal_a] < alpha)       # True causal called significant
# FP_a <- sum(combined_results$p[!combined_results$causal_a] < alpha)      # Non-causal called significant
# FN_a <- sum(combined_results$p[combined_results$causal_a] >= alpha)      # True causal NOT called significant
# TN_a <- sum(combined_results$p[!combined_results$causal_a] >= alpha)     # Non-causal NOT called significant
# 
# 
# 
# TP_d <- sum(combined_results$p[combined_results$causal_d] < alpha)       # True causal called significant
# FP_d <- sum(combined_results$p[!combined_results$causal_d] < alpha)      # Non-causal called significant
# FN_d <- sum(combined_results$p[combined_results$causal_d] >= alpha)      # True causal NOT called significant
# TN_d <- sum(combined_results$p[!combined_results$causal_d] >= alpha)     # Non-causal NOT called significant
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

# sum(combined_results$P < alpha)
