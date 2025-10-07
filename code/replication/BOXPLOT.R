# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)

# Function to process one scenario
process_scenario <- function(results_dir, pheno_dir, scenario_label) {
  p_folders <- c("10E-3", "10E-5", "10E-8")
  ld_thresholds <- c(0.1, 0.3, 0.5)
  n_reps <- 20
  
  # Lists to store correlation and regression results
  all_results_corr <- list()
  all_results_reg <- list()
  
  for (p_folder in p_folders) {
    p_dir <- file.path(results_dir, p_folder)
    if (!dir.exists(p_dir)) stop("Directory not found: ", p_dir)
    
    corr_list <- list()
    reg_list <- list()  # To store regression results
    
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
      
      merged_scan <- merge(pheno, prs_scan, by = "IID")
      merged_gwas <- merge(pheno, prs_gwas, by = "IID")
      
      # === Correlation Calculations (existing) ===
      scan_corrs <- sapply(ld_thresholds, function(ld) {
        col <- paste0("r2_", ld)
        cor(merged_scan$pheno, merged_scan[[col]], use = "complete.obs")
      })
      gwas_corrs <- sapply(ld_thresholds, function(ld) {
        col <- paste0("r2_", ld)
        cor(merged_gwas$pheno, merged_gwas[[col]], use = "complete.obs")
      })
      
      res_rep_corr <- data.table(
        rep = i,
        ld_threshold = factor(ld_thresholds, levels = c(0.1, 0.3, 0.5)),
        prs_scan_r = scan_corrs,
        prs_gwas_r = gwas_corrs
      )
      corr_list[[i]] <- res_rep_corr
      
      
      # === Regression: pheno ~ prs (extract intercept and slope) ===
      reg_rep_list <- list()
      for (ld in ld_thresholds) {
        col <- paste0("r2_", ld)
        
        fit_scan <- lm(as.formula(paste("pheno ~", col)), data = merged_scan)
        coef_scan <- coef(fit_scan)
        
        # Regression for gwas
        fit_gwas <- lm(as.formula(paste("pheno ~", col)), data = merged_gwas)
        coef_gwas <- coef(fit_gwas)
        
        reg_rep <- data.table(
          rep = i,
          ld_threshold = factor(ld, levels = c(0.1, 0.3, 0.5)),
          method = "scan",
          intercept = coef_scan[1],
          slope = coef_scan[2]
        )
        
        reg_rep_gwas <- data.table(
          rep = i,
          ld_threshold = factor(ld, levels = c(0.1, 0.3, 0.5)),
          method = "gwas",
          intercept = coef_gwas[1],
          slope = coef_gwas[2]
        )
        
        reg_rep_list[[ld]] <- rbind(reg_rep, reg_rep_gwas)
      }
      
      # Combine all LD thresholds for this rep
      reg_list[[i]] <- rbindlist(reg_rep_list)
    }
    
    # Combine all reps for this p_threshold
    res_all_reps_corr <- rbindlist(corr_list)
    res_all_reps_reg  <- rbindlist(reg_list)
    
    # Add p_threshold column
    res_all_reps_corr[, p_threshold := p_folder]
    res_all_reps_reg[,  p_threshold := p_folder]
    
    all_results_corr[[p_folder]] <- res_all_reps_corr
    all_results_reg[[p_folder]]  <- res_all_reps_reg
  }
  
  # Combine all p-thresholds
  combined_corr <- rbindlist(all_results_corr)
  combined_reg  <- rbindlist(all_results_reg)
  
  # Compute R² from r
  combined_corr[, `:=`(
    prs_scan_r2 = prs_scan_r^2,
    prs_gwas_r2 = prs_gwas_r^2
  )]
  
  # Add scenario label
  combined_corr[, Scenario := scenario_label]
  combined_reg[,  Scenario := scenario_label]
  
  # Optional: Set key for clarity
  setkeyv(combined_reg, c("rep", "ld_threshold", "p_threshold", "method", "Scenario"))
  
  # Return both results as a list
  return(list(
    correlations = combined_corr,
    regressions  = combined_reg
  ))
}
# --- Load All Three Scenarios ---
scenario1 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario1/pheno",
  scenario_label = "Scenario 1: Pure Additive"
)

scenario2 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario2/pheno",
  scenario_label = "Scenario 2: Partial Additive with incomplete dominance"
)

scenario3 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario3/pheno",
  scenario_label = "Scenario 3: Partial Additive with complete dominance"
)

scenario4 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno",
  scenario_label = "Scenario 4: Partial Additive with complete dominance"
)
scenario5 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario5",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario5/pheno",
  scenario_label = "Scenario 5: Partial Additive with complete dominance"
)
scenario6 <- process_scenario(
  results_dir = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario6",
  pheno_dir   = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario6/pheno",
  scenario_label = "Scenario 6: Partial Additive with complete dominance"
)
# --- Load and extract correlations only ---
scenario1_corr <-scenario1$correlations

# scenario2_corr <- scenario2$correlations

scenario3_corr <- scenario3$correlations

scenario4_corr <- scenario4$correlations

scenario5_corr <- scenario5$correlations
scenario6_corr <- scenario6$correlations
# Combine all three scenarios
plot_data <- rbindlist(list(scenario4_corr, scenario5_corr), fill = TRUE)

# Reshape to long format for plotting
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

# Create composite x-axis label: P-value + LD threshold
plot_data_long <- plot_data_long %>%
  mutate(
    x_label = paste0("P: ", p_threshold, "\nLD: ", ld_threshold)
  ) %>%
  # Ensure correct order: sort by p_threshold, then ld_threshold
  mutate(
    x_label = factor(x_label, levels = unique(x_label[order(p_threshold, as.numeric(as.character(ld_threshold)))]))
  )

p <- ggplot(plot_data_long, aes(x = x_label, y = Correlation, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ Scenario, nrow = 4, labeller = label_wrap_gen(width = 60)) +  # nrow = 4 now
  labs(
    title = "PRS Performance Across Four Simulation Scenarios",
    subtitle = "Comparison of PRS-Scan vs PRS-GWAS in Predicting Phenotype (R²)",
    x = "P-value Threshold and LD Clumping Threshold",
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

# Save the static ggplot (before making it interactive)
ggsave(
  filename = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/prs_comparison_plot_3Scenarios.png",  # or .pdf, .jpg, etc.
  plot = p,                              # the ggplot object
  device = "png",                        # optional: specify device
  width = 12,                            # adjust as needed
  height = 8,
  dpi = 300,
  units = "in"
)

# Convert to interactive plot
p_interactive <- ggplotly(p, tooltip = c("x", "y", "fill", "Scenario"))

# Display interactive plot
p_interactive



## Regression
# --- Extract regression results for each scenario ---
scenario1_reg <- scenario1$regressions
scenario2_reg <- scenario2$regressions
scenario3_reg <- scenario3$regressions
scenario4_reg <- scenario4$regressions

# Combine all regression results
reg_plot_data <- rbindlist(
  list(scenario1_reg, scenario2_reg),
  fill = TRUE
)

# Reshape to long format
reg_long <- melt(reg_plot_data,
                 measure.vars = c("intercept", "slope"),
                 variable.name = "coefficient",
                 value.name = "value"
)

# Fix labels
reg_long[, method_label := ifelse(method == "scan", "PRS-SCAN", "PRS-GWAS")]
# Optional: shorten scenario names if too long
reg_plot_data[, Scenario_short := gsub("Scenario (\\d):.*", "Scenario \\1", Scenario)]

# Create method labels
reg_plot_data[, method_label := ifelse(method == "scan", "PRS-SCAN", "PRS-GWAS")]
reg_plot_data[, method_label := factor(method_label, levels = c("PRS-SCAN", "PRS-GWAS"))]
reg_plot_data[, Scenario_short := factor(Scenario_short, levels = unique(Scenario_short))]

# Define colors
method_colors <- c("PRS-SCAN" = "red", "PRS-GWAS" = "blue")

# Unique LD thresholds
ld_levels <- unique(reg_plot_data$ld_threshold)

# List to store plots
plot_list <- list()

# Loop over each LD threshold
for (ld in ld_levels) {
  # Subset data for this LD threshold
  dat_ld <- reg_plot_data[ld_threshold == ld]
  
  # Reshape to long format for slope and intercept
  dat_long <- melt(
    dat_ld,
    id.vars = c("rep", "method_label", "Scenario", "Scenario_short"),
    measure.vars = c("slope", "intercept"),
    variable.name = "coefficient",
    value.name = "value"
  )
  
  # Format coefficient labels
  dat_long[, coefficient := ifelse(coefficient == "slope", "Slope", "Intercept")]
  
  # Create the plot
  p <- ggplot(dat_long, aes(x = Scenario_short, y = value, fill = method_label)) +
    geom_boxplot(alpha = 0.8, width = 0.6, outlier.alpha = 0.4) +
    scale_fill_manual("PRS Method", values = method_colors) +
    facet_wrap(~ coefficient, scales = "free_y", ncol = 1) +
    labs(
      title = paste("Regression Coefficients — LD Threshold: r² <", ld),
      x = "Simulation Scenario",
      y = "Coefficient Value"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) +
    # Add reference lines
    geom_hline(data = data.table(coefficient = c("Slope", "Intercept"),
                                 hline_val = c(1, 0)),
               aes(yintercept = hline_val), linetype = "dashed", color = "gray50", size = 0.8)
  
  # Store plot
  plot_list[[as.character(ld)]] <- p
}

# Print all plots
plot_list
