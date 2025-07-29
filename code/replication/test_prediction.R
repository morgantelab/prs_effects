# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
rep <- as.integer(args[1])
results_dir <- args[2]
ld_thresholds <- args[3]
test_ID <- args[4]
effects_file <- args[5]
clump_file <- args[6]
geno_pattern <- args[6]
plink_path <- args[7]
output_file <- args[8]

cat("\nStarting PRS calculation for replicate", rep, "\n")
cat("Results directory:", results_dir, "\n")
cat("LD thresholds:", ld_thresholds, "\n")
cat("TEST ID files", test_ID, "\n")
cat("Effect Files:", effects_file, "\n")
cat("Clump Files:", clump_file, "\n")
cat("Geno File Names:", geno_pattern, "\n")
cat("PLINK path", plink_path, "\n")

library(data.table)
library(genio)

ld_thresholds <- as.numeric(unlist(strsplit(ld_thresholds, ",")))
cat("LD thresholds", ld_thresholds, "\n")

# Load test set IDs
test_ids <- fread(test_ID, header = FALSE)[[1]]

# Load effect estimates
eff <- fread(effects_file)

# Prepare PRS storage
prs_matrix <- matrix(0, nrow = length(test_ids), ncol = length(ld_thresholds))
colnames(prs_matrix) <- sprintf("r2_%s", format(ld_thresholds, nsmall = 1))
temp_dir <- file.path(results_dir, sprintf("temp_rep%d", rep))

for (r2_idx in seq_along(ld_thresholds)) {
  r2 <- ld_thresholds[r2_idx]
  r2_str <- format(r2, nsmall = 1)
  cat("  LD threshold r² =", r2_str, "\n")
  
  # Load clumped SNPs
  clumped <- fread(clump_file)
  setnames(clumped, c("ID", "#CHROM"), c("SNP", "chromosome"))
  
  # Get effect sizes for clumped SNPs
  clumped_effects <- merge(clumped[, .(SNP, chromosome)], 
                           eff, 
                           by = c("SNP", "chromosome"))
  setnames(clumped_effects, "chromosome", "CHR")


  for (chr in 1:22) {
    cat("    Processing chromosome", chr, "\n")
    
    # Get SNPs for this chromosome
    chr_snps <- clumped_effects[CHR == chr, 
                                .(SNP, ref_allele, alt_allele, a_hat, d_hat)]
    snp_list_file <- file.path(temp_dir, sprintf("chr%d_snps.txt", chr))
    fwrite(chr_snps[, .(SNP)], snp_list_file, col.names = FALSE)  
    # Prepare genotype files
    plink_file <- paste0(geno_pattern, chr, "_filtered" )
    system(sprintf(
      "%s --bfile %s --keep %s --extract %s --make-bed --out %s --memory 2000 --threads 10",
      plink_path,
      plink_prefix,
      test_ids,
      snp_list_file,
      temp_prefix
    ))

    # Read genotype data
    plink_data <- genio::read_plink(temp_prefix, verbose = FALSE)
    fam_ids <- plink_data$fam$id
    snp_ids <- plink_data$bim$id
    X <- plink_data$X  
    matched_snp_ids <- snp_ids[snp_logical]

    # Match SNP effects to the order in the genotype matrix
    snp_effects <- chr_snps[match(matched_snp_ids, chr_snps$SNP), ]
    
    # Calculate PRS contributions
    for (j in 1:ncol(X)) {
      g <- X[, j]
      a <- snp_effects$a_hat[j]
      d <- snp_effects$d_hat[j]
      
      # Direct genotype-based calculation
      contrib <- rep(0, length(g))
      contrib[g == 0] <- -a  # Homozygous reference
      contrib[g == 1] <- d   # Heterozygous
      contrib[g == 2] <- a   # Homozygous alternate
      
      # Add to PRS matrix
      prs_matrix[, r2_idx] <- prs_matrix[, r2_idx] + contrib
    }
  }
}
# Cleanup temporary files
unlink(temp_dir, recursive = TRUE)

# Save PRS results
prs_df <- data.table(iid = test_ids, prs_matrix)
fwrite(prs_df, output_file, sep = "\t")
cat("Saved PRS results to:", output_file, "\n")