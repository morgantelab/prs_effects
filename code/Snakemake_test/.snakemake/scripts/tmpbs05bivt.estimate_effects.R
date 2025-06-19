
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list(),
    output = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/results/chr_effects_test/chr5_effects.txt', "chr_effects" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results/chr_effects_test/chr5_effects.txt'),
    params = list('5', '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter', '/data2/morgante_lab/ukbiobank_projects/prs_effects/results', "chr" = '5', "base_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter', "results_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results'),
    wildcards = list('5', "chr" = '5'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("base_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter', "results_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results', "plink2_path" = '/data2/morgante_lab/yuxuans/software/plink2', "log_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/logs', "error_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/errors', "code_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/code', "R_LIB1" = '/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.2', "R_LIB2" = '/opt/ohpc/pub/Software/R-4.2.3/lib64/R/library', "cluster_config" = list("process_chromosome" = list("mem_mb" = 768000, "time" = '3-00:00:00', "partition" = 'fm-bigmem-2', "threads" = 4, "mail_user" = 'yuxuans@clemson.edu', "mail_type" = 'FAIL'), "combine_results" = list("mem_mb" = 51200, "time" = '1-00:00:00', "partition" = 'fm-bigmem-2', "threads" = 4, "mail_user" = 'yuxuans@clemson.edu', "mail_type" = 'FAIL')), "cluster" = 'mkdir -p {log_dir}/{rule} && sbatch --partition={resources.partition} --cpus-per-task={threads} --mem={resources.mem_mb} --job-name={rule}_{wildcards} --output={log_dir}/{rule}/%j_{rule}_{wildcards}.out --error={error_dir}/{rule}/%j_{rule}_{wildcards}.err --time={resources.time} --mail-user={resources.mail_user} --mail-type={resources.mail_type}'),
    rule = 'process_chromosome',
    bench_iteration = as.numeric(NA),
    scriptdir = '/data2/morgante_lab/ukbiobank_projects/prs_effects/code',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(genio)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
chr <- as.integer(args[1])
base_dir <- args[2]
results_dir <- args[3]

# Create output directory
chr_dir <- file.path(results_dir, "chr_effects_test")
if (!dir.exists(chr_dir)) dir.create(chr_dir, recursive = TRUE)

print(.libPaths())
###----------------------- Effect Calculation Function ------------------------
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
        (var_aa/(4*n_aa))^2/(n_aa-1))
    t_d <- d_hat/sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}

###----------------------- Process Single Chromosome --------------------------
cat("Processing chromosome", chr, "\n")

# Read phenotype data
pheno_file <- file.path(results_dir, "simulated_pheno_large.txt")
pheno <- fread(pheno_file,
               colClasses = c("FID" = "character", "IID" = "character", "pheno" = "numeric"))

# Read PLINK files
plink_file <- file.path(base_dir, paste0("ukb_chr", chr, "_filtered"))
plink_data <- read_plink(plink_file)

# Create genotype data frame
geno_df <- as.data.frame(t(plink_data$X))
colnames(geno_df) <- plink_data$bim$id
geno_df <- cbind(plink_data$fam[, c("fam", "id")], geno_df)
names(geno_df)[1:2] <- c("FID", "IID")

# Extract SNP info
snp_info <- plink_data$bim
names(snp_info) <- c("chromosome", "SNP", "pos_cM", "position", "ref_allele", "alt_allele")
snp_info <- snp_info[, c("chromosome", "SNP", "position", "ref_allele", "alt_allele")]

# Merge with phenotype data
full_data <- merge(pheno, geno_df, by = c("FID", "IID"))
snp_cols <- setdiff(colnames(full_data), c("FID", "IID", "pheno"))

# Calculate Welch statistics
chr_results <- calculate_welch(full_data, snp_cols)
chr_results$chr <- chr

# Add SNP information
chr_results <- merge(chr_results, snp_info, by = "SNP", all.x = TRUE) %>%
  dplyr::select(chromosome, position, SNP, ref_allele, alt_allele, everything())

# Save results
output_file <- file.path(chr_dir, paste0("chr", chr, "_effects.txt"))
fwrite(chr_results, output_file, sep = "\t", quote = FALSE)
cat("Saved results for chromosome", chr, "\n")