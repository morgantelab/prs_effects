# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)

# Function to process one scenario
process_scenario_simple <- function(results_dir, pheno_dir, scenario_label) {
  n_reps <- 20
  
  all_results <- list()
  
  for (i in 1:n_reps) {
    pheno_file <- file.path(pheno_dir, paste0("simulated_pheno_rep", i, ".txt"))
    scan_file  <- file.path(results_dir, paste0("test_prs_scan_rep", i, ".txt"))
    gwas_file  <- file.path(results_dir, paste0("test_prs_gwas_rep", i, ".txt"))
    
    true_pheno    <- fread(pheno_file)
    prs_scan <- fread(scan_file)
    prs_gwas <- fread(gwas_file)
    
    merged_scan <- merge(true_pheno, prs_scan, by = "IID")
    merged_gwas <- merge(true_pheno, prs_gwas, by = "IID")
    
    # Compute correlation
    scan_corr <- cor(merged_scan$pheno, merged_scan$PRS)^2
    gwas_corr <- cor(merged_gwas$pheno, merged_gwas$PRS_total)^2
    
    res <- data.table(
      replicate = i,
      prs_method = c("PRS-Scan", "PRS-GWAS"),
      R2 = c(scan_corr^2, gwas_corr^2),
      Scenario = scenario_label
    )
    
    all_results[[i]] <- res
  }
  
  combined <- rbindlist(all_results)
  return(combined)
}

# --- Process scenarios ---
scenario1 <- process_scenario_simple(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1/pheno",
  scenario_label = "Scenario 1"
)

scenario2 <- process_scenario_simple(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2/pheno",
  scenario_label = "Scenario 2"
)

scenario3 <- process_scenario_simple(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3/pheno",
  scenario_label = "Scenario 3"
)

scenario4 <- process_scenario_simple(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno",
  scenario_label = "Scenario 4"
)


# Combine all scenarios
plot_data <- rbindlist(list(scenario1, scenario2, scenario3, scenario4), fill = TRUE)

# --- Plot ---
p <- ggplot(plot_data, aes(x = prs_method, y = R2, fill = prs_method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Scenario, nrow = 4) +
  labs(
    title = "PRS Performance Across Simulation Scenarios",
    subtitle = "R² of PRS-Scan vs PRS-GWAS with Phenotype",
    x = "PRS Method",
    y = "R² with Phenotype",
    fill = "PRS Method"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "gray80", size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 12, color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))

p

ggsave(
  filename = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/PRS_performance_causal.png",
  plot = p,
  width = 8,
  height = 6
)
