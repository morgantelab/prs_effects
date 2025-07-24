#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

# Set paths
results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication"
ldetect_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR"

# Process each replicate
for (rep in 1:20) {
  cat("\nProcessing replicate", rep, "\n")
  
  # Input and output files
  effects_file <- file.path(results_dir, sprintf("effects_rep%d.txt", rep))
  output_file <- file.path(results_dir, sprintf("selected_snps_rep%d.txt", rep))
  
  # Read effects data
  if (!file.exists(effects_file)) {
    cat("Effects file not found:", effects_file, "\n")
    next
  }
  
  effects <- fread(effects_file)
  
  # Calculate replicate-specific Bonferroni threshold
  total_snps <- nrow(effects)
  bonferroni_threshold <- 0.05 / total_snps
  cat("Total SNPs in replicate", rep, ":", total_snps, "\n")
  cat("Bonferroni threshold:", bonferroni_threshold, "\n")
  
  selected_snps <- data.table()
  
  # Process each chromosome
  for (chr in 1:22) {
    cat("Processing chromosome", chr, "\n")
    
    # Read LD block boundaries
    bed_file <- file.path(ldetect_dir, sprintf("fourier_ls-chr%d.bed", chr))
    if (!file.exists(bed_file)) {
      cat("LDetect file not found:", bed_file, "\n")
      next
    }
    
    ldetect <- fread(bed_file, col.names = c("chr", "start", "end", "name"))
    
    # Filter effects for current chromosome
    chr_effects <- effects[chromosome == chr]
    if (nrow(chr_effects) == 0) next
    
    # Convert SNP positions to 0-based coordinates
    chr_effects[, pos0 := position - 1]
    
    # Assign SNPs to LD blocks
    block_index <- findInterval(chr_effects$pos0, ldetect$start)
    final_block_index <- rep(NA_integer_, nrow(chr_effects))
    
    for (i in 1:nrow(ldetect)) {
      idx <- which(block_index == i)
      if (length(idx) == 0) next
      in_block <- chr_effects$pos0[idx] < ldetect$end[i]
      final_block_index[idx[in_block]] <- i
    }
    
    chr_effects[, block_idx := final_block_index]
    
    # Process each block
    for (block_id in 1:nrow(ldetect)) {
      block_snps <- chr_effects[block_idx == block_id]
      if (nrow(block_snps) == 0) next
      
      # Filter significant SNPs (using additive p-value)
      sig_snps <- block_snps[p_a < bonferroni_threshold]
      if (nrow(sig_snps) == 0) next
      
      # Select SNP with the highest p-value (least significant among significant)
      selected <- sig_snps[which.max(p_a)]
      selected_snps <- rbind(selected_snps, selected)
    }
  }
  
  # Save results
  if (nrow(selected_snps) > 0) {
    # Remove temporary columns
    selected_snps[, c("pos0", "block_idx") := NULL]
    fwrite(selected_snps, output_file, sep = "\t")
    cat("Saved", nrow(selected_snps), "selected SNPs to", output_file, "\n")
  } else {
    cat("No significant SNPs found for replicate", rep, "\n")
  }
}

cat("\nAll replicates processed.\n")