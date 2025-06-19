
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
    input = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/results/simulated_pheno_large.txt', "pheno" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results/simulated_pheno_large.txt'),
    output = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/temp/chr20_effects_test.txt'),
    params = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr20_filtered', "plink_base" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr20_filtered'),
    wildcards = list('20', "chr" = '20'),
    threads = 1,
    log = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/process_chromosome/chr20.out', '/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/process_chromosome/chr20.err', "out" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/process_chromosome/chr20.out', "err" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/process_chromosome/chr20.err'),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', "mem_mb" = 1000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp'),
    config = list("base_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter', "results_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results', "temp_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/temp', "log_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/logs', "code_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/code', "error_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/errors', "plink2_path" = '/data2/morgante_lab/yuxuans/software/plink2', "__default__" = list("account" = 'fm-bigmem-2', "partition" = 'fm-bigmem-2', "mail-user" = 'yuxuans@clemson.edu', "mail-type" = 'ALL', "mem" = '800G', "time" = '12:00:00'), "process_chromosome" = list("mem" = '800G', "time" = '12:00:00'), "combine_results" = list("mem" = '200G', "time" = '01:00:00'), "cluster" = 'mkdir -p {log_dir}/{rule} && sbatch --account={resources.account} --partition={resources.partition} --cpus-per-task={threads} --mem={resources.mem} --job-name={rule}_{wildcards} --output={log_dir}/{rule}/%j_{rule}_{wildcards}.out --error={error_dir}/{rule}/%j_{rule}_{wildcards}.err --time={resources.time} --mail-user={resources.mail-user} --mail-type={resources.mail-type}', "latency-wait" = 60, "printshellcmds" = TRUE),
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

# Get Snakemake input and output paths
pheno_file <- snakemake@input[["pheno"]]
plink_base <- snakemake@params[["plink_base"]]
output_file <- snakemake@output[[1]]
chr <- snakemake@wildcards[["chr"]]
# Set fixed paths
base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects"
results_dir <- file.path(base_dir, "results")
temp_dir <- file.path(base_dir, "temp")
log_dir <- file.path(base_dir, "logs")

# Load required packages

library(data.table)
library(dplyr)
library(genio)

# Read phenotype data
cat("Read in Phenotype files", "\n")

pheno <- fread(pheno_file, colClasses = c("character", "character", "numeric"))

# Welch's t-test function
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
        (var_aa/(4*n_aa))^2/(n_aa-1)
    )
    t_d <- d_hat/sqrt(var_d)
    p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)
    
    results[j,] <- list(snp, a_hat, d_hat, p_a, p_d)
  }
  
  return(results)
}

# Process chromosome data
cat("PLINK base:", plink_base, "\n")

plink_data <- read_plink(plink_base)
cat("Successfully read in plink data", "\n")

geno_df <- as.data.frame(t(plink_data$X))
cat("Successfully transpose Geno", "\n")
colnames(geno_df) <- plink_data$bim$id
geno_df <- cbind(plink_data$fam[c("fam", "id")], geno_df)
names(geno_df)[1:2] <- c("FID", "IID")
snp_info <- plink_data$bim %>%
  dplyr::select(id, chr, pos, ref, alt) %>%
  dplyr::rename(SNP = id, 
                chromosome = chr,
                position = pos,
                ref_allele = ref,
                alt_allele = alt)
cat("Successfully Extract Geno Information", "\n")

# Merge phenotype and genotype data
full_data <- merge(pheno, geno_df, by = c("FID", "IID"))
snp_cols <- setdiff(colnames(full_data), c("FID", "IID", "pheno"))

# Calculate effects
chr_results <- calculate_welch(full_data, snp_cols)
chr_results$chr <- chr

# Add SNP information
chr_results <- chr_results %>%
  left_join(snp_info, by = "SNP") %>%
  dplyr::select(chromosome, position, SNP, ref_allele, alt_allele, everything())

# Write output
output_file <- file.path(temp_dir, paste0("chr", chr, "_effects_test.txt"))
cat("Writing results to:", output_file, "\n")

write.table(chr_results, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Completed processing chromosome", chr, "\n")