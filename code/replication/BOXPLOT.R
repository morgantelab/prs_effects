library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Function to process one scenario
process_scenario <- function(results_dir, pheno_dir, scenario_label) {
  p_folders <- c("10E-3", "10E-5", "10E-8")
  ld_thresholds <- c(0.1, 0.3, 0.5)
  n_reps <- 20
  all_results <- list()
  
  for (p_folder in p_folders) {
    p_dir <- file.path(results_dir, p_folder)
    if (!dir.exists(p_dir)) stop("Directory not found: ", p_dir)
    
    res_list <- list()
    
    for (i in 1:n_reps) {
      pheno_file <- file.path(pheno_dir, paste0("simulated_pheno_rep", i, ".txt"))
      scan_file  <- file.path(p_dir, paste0("test_prs_scan_rep", i, ".txt"))
      gwas_file  <- file.path(p_dir, paste0("test_prs_gwas_rep", i, ".txt"))
      
      if (!file.exists(pheno_file)) stop("Not found: ", pheno_file)
      if (!file.exists(scan_file))  stop("Not found: ", scan_file)
      if (!file.exists(gwas_file))  stop("Not found: ", gwas_file)
      
      pheno    <- fread(pheno_file, header = TRUE)
      prs_scan <- fread(scan_file, header = TRUE)
      prs_gwas <- fread(gwas_file, header = TRUE)
      
      merged_scan <- merge(pheno, prs_scan, by.x = "IID", by.y = "iid")
      merged_gwas <- merge(pheno, prs_gwas, by = "IID")
      
      scan_corrs <- sapply(ld_thresholds, function(ld) {
        col <- paste0("r2_", ld)
        cor(merged_scan$pheno, merged_scan[[col]], use = "complete.obs")
      })
      gwas_corrs <- sapply(ld_thresholds, function(ld) {
        col <- paste0("r2_", ld)
        cor(merged_gwas$pheno, merged_gwas[[col]], use = "complete.obs")
      })
      
      res_rep <- data.table(
        rep = i,
        ld_threshold = factor(ld_thresholds, levels = c(0.1, 0.3, 0.5)),
        prs_scan_r = scan_corrs,
        prs_gwas_r = gwas_corrs
      )
      res_list[[i]] <- res_rep
    }
    
    # Combine all reps for this p_threshold
    res_all_reps <- rbindlist(res_list)
    
    # ✅ Add p_threshold column (only once!)
    res_all_reps[, p_threshold := p_folder]
    
    all_results[[p_folder]] <- res_all_reps
  }
  
  # Combine all p-thresholds
  combined <- rbindlist(all_results)  # No idcol!
  
  # Compute R²
  combined[, `:=`(
    prs_scan_r2 = prs_scan_r^2,
    prs_gwas_r2 = prs_gwas_r^2
  )]
  
  # Add scenario label
  combined[, Scenario := scenario_label]
  
  return(combined)
}
# --- Load Both Scenarios ---
scenario1 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1/pheno",
  scenario_label = "Scenario 1: Pure Additive"
)

scenario2 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2/pheno",
  scenario_label = "Scenario 2: Partial Additive"
)

# Combine using data.table's rbindlist (safe and clean)
plot_data <- rbindlist(list(scenario1, scenario2), fill = TRUE)

# Now proceed with dplyr
plot_data_long <- plot_data %>%
  select(Scenario, p_threshold, ld_threshold, rep, prs_scan_r2, prs_gwas_r2) %>%
  pivot_longer(
    cols = c(prs_scan_r2, prs_gwas_r2),
    names_to = "Method",
    values_to = "R2"
  ) %>%
  mutate(
    Method = case_when(
      Method == "prs_scan_r2" ~ "PRS-Scan",
      Method == "prs_gwas_r2" ~ "PRS-GWAS"
    ),
    Method = factor(Method, levels = c("PRS-Scan", "PRS-GWAS")),
    p_threshold = factor(p_threshold, levels = c("10E-3", "10E-5", "10E-8"))
  ) %>%
  rename(Correlation = R2)

# Create x-axis label: P-value + LD threshold
plot_data_long <- plot_data_long %>%
  mutate(
    x_label = paste0("P: ", p_threshold, "\nLD: ", ld_threshold)
  ) %>%
  # Keep order: by P-threshold, then LD
  mutate(
    x_label = factor(x_label, levels = unique(x_label[order(p_threshold, ld_threshold)]))
  )

ggplot(plot_data_long, aes(x = x_label, y = Correlation, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Scenario, nrow = 2, labeller = label_wrap_gen(width = 50)) +
  labs(
    title = "PRS Performance: Scenario 1 vs Scenario 2",
    x = "P-value Threshold and LD Clumping Threshold",
    y = "R² with Phenotype",
    fill = "Method"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))
