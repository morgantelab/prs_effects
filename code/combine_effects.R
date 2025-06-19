#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[1]

# Set directory paths
chr_dir <- file.path(results_dir, "chr_effect")
# List all chromosome files
chr_files <- list.files(chr_dir, pattern = "chr\\d+_effects\\.txt", full.names = TRUE)
cat("Found", length(chr_files), "chromosome files\n")

# Read and combine
combined <- rbindlist(lapply(chr_files, fread), fill = TRUE)

# Save combined results
output_file <- file.path(results_dir, "effect_size_large_test.txt")
fwrite(combined, output_file, sep = "\t", quote = FALSE)
cat("Combined results saved to:", output_file, "\n")