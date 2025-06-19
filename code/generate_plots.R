#!/usr/bin/env Rscript

# Set fixed paths
base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects"
results_dir <- file.path(base_dir, "results")

# Input and output files
input_file <- file.path(results_dir, "effect_size_large_test.txt")
output_add <- file.path(results_dir, "manhattan_additive_large_test.png")
output_dom <- file.path(results_dir, "manhattan_dominance_large_test.png")

effect_results <- fread(input_file)

# Manhattan plot function
create_manhattan <- function(data, p_col, title) {
  # Prepare data
  manhattan_data <- data %>%
    rename(P = !!p_col) %>%
    rename(CHR = chromosome) %>%
    rename(BP = position) %>%
    filter(!is.na(P))
  
  # Compute cumulative positions
  data_cum <- manhattan_data %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(BP)) %>% 
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(manhattan_data, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot)
  
  # Calculate axis positions
  axisdf <- data_cum %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  # Bonferroni correction
  bonferroni <- 0.05 / nrow(data_cum)
  
  # Create plot
  ggplot(data_cum, aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR %% 2)), alpha = 0.7, size = 0.8) +
    scale_color_manual(values = c("gray60", "gray30")) +
    scale_x_continuous(labels = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_hline(yintercept = -log10(bonferroni), color = "red", linetype = "dashed") +
    labs(
      title = title,
      x = "Chromosome",
      y = expression(-log[10](p-value))
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# Generate and save plots
cat("Generating Manhattan plots...\n")
cat("Additive effects plot:", output_add, "\n")
cat("Dominance effects plot:", output_dom, "\n")

tryCatch({
  p_add <- create_manhattan(effect_results, "p_a", "Additive Effects")
  ggsave(output_add, p_add, width = 12, height = 6, dpi = 300)
  
  p_dom <- create_manhattan(effect_results, "p_d", "Dominance Effects")
  ggsave(output_dom, p_dom, width = 12, height = 6, dpi = 300)
  
  cat("Successfully generated Manhattan plots\n")
}, error = function(e) {
  stop("Error generating plots: ", e$message)
})