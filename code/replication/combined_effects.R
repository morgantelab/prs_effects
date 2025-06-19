#!/usr/bin/env Rscript
# combine_effects.R: Combine chromosome effect files into a single output

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chr_dir <- args[1]
output_file <- args[2]
rep <- args[3]

# Load required packages
library(data.table)
library(dplyr)

# Find all chromosome result files with consistent pattern
input_files <- list.files(chr_dir, 
                          pattern = sprintf("^chr\\d{1,2}/rep%s.txt$", rep), 
                          full.names = TRUE)

# Read and combine files efficiently
cat("Combining", length(input_files), "chromosome files...\n")
combined_data <- rbindlist(lapply(input_files, function(f) {
  cat("Processing:", f, "\n")
  tryCatch({
    dt <- fread(f)
    # Add basic validation
    if (nrow(dt) == 0) {
      warning("Empty file: ", f)
      return(NULL)
    }
    required_cols <- c("chromosome", "position", "SNP", "ref_allele", 
                       "alt_allele", "a_hat", "d_hat", "p_a", "p_d")
    if (!all(required_cols %in% names(dt))) {
      warning("Missing columns in: ", f)
      return(NULL)
    }
    return(dt)
  }, error = function(e) {
    warning("Error reading ", f, ": ", e$message)
    return(NULL)
  })
}), use.names = TRUE, fill = TRUE)

# Final data checks
if (nrow(combined_data) == 0) {
  stop("No data combined! Check input files.")
}

# Write output
cat("Writing combined output to:", output_file, "\n")
fwrite(combined_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
