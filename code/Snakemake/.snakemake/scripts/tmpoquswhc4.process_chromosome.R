
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
    input = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.bed', '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.bim', '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.fam', '/data2/morgante_lab/ukbiobank_projects/prs_effects/results/simulated_pheno_large.txt', "bed" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.bed', "bim" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.bim', "fam" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered.fam', "pheno" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results/simulated_pheno_large.txt'),
    output = list('/data2/morgante_lab/ukbiobank_projects/prs_effects/temp/chr22_effects_test.txt'),
    params = list(),
    wildcards = list('22', "chr" = '22'),
    threads = 1,
    log = list(),
    resources = list('mem_mb', 'mem_mib', 'disk_mb', 'disk_mib', 'tmpdir', 'time_min', 'partition', 'mem', "mem_mb" = 800000, "mem_mib" = 954, "disk_mb" = 1000, "disk_mib" = 954, "tmpdir" = '/tmp', "time_min" = '12:00:00', "partition" = 'fm-bigmem-2', "mem" = '800G'),
    config = list("base_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/filter', "results_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/results', "temp_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/temp', "log_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/logs', "code_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/code', "error_dir" = '/data2/morgante_lab/ukbiobank_projects/prs_effects/errors', "plink2_path" = '/data2/morgante_lab/yuxuans/software/plink2', "resources" = c('cpus=400', 'mem_mb=8000000', 'jobweight=40', 'solo=1'), "cluster" = 'mkdir -p {log_dir}/{rule} && sbatch --account={resources.account} --partition={resources.partition} --ntasks-per-core={resources.ntpc} --ntasks={resources.tasks} --cpus-per-task={threads} --mem={resources.mem} --job-name={rule}_{wildcards} --output={log_dir}/{rule}/%j_{rule}_{wildcards}.out --error={error_dir}/{rule}/%j_{rule}_{wildcards}.err --time={resources.time} --mail-user={resources.mail-user} --mail-type={resources.mail-type}', "default-resources" = c("partition='fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4'", 'mem_mb_per_cpu=12000', 'time="4-00:00:00"', 'nodes=4', 'tasks=4', 'ntpc=4', 'mem_mb=10000'), "latency-wait" = 60, "printshellcmds" = TRUE, "__default__" = list("account" = 'fm-bigmem-2', "partition" = 'fm-bigmem-2', "output" = '{log_dir}/Snakemake_{rule}_{wildcards}.%j.out', "error" = '{error_dir}/Snakemake_{rule}_{wildcards}.%j.err', "mail-user" = 'yuxuans@clemson.edu', "mail-type" = 'ALL', "mem" = '8G', "time" = '01:00:00'), "process_chromosome" = list("account" = 'fm-bigmem-2', "partition" = 'fm-bigmem-2', "mem" = '800G', "time" = '12:00:00', "cpus-per-task" = 4), "combine_results" = list("account" = 'fm-bigmem-2', "partition" = 'fm-bigmem-2', "mem" = '200G', "time" = '01:00:00'), "generate_plots" = list("account" = 'fm-bigmem-2', "partition" = 'fm-bigmem-2', "mem" = '300G', "time" = '02:00:00')),
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
args <- commandArgs(trailingOnly = TRUE)

chr <- args[1]
bed_file <- args[2]
bim_file <- args[3]
fam_file <- args[4]

# Set fixed paths
base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects"
results_dir <- file.path(base_dir, "results")
temp_dir <- file.path(base_dir, "temp")
log_dir <- file.path(base_dir, "logs")
code_dir <- file.path(base_dir, "code/Snakemake")

.libPaths(c(
  "/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.1",
  "/opt/ohpc/pub/libs/gnu9/R/4.1.2/lib64/R/library"
))

# Construct PLINK base name by removing .bed extension
plink_base <- sub("\\.bed$", "", bed_file)

# Load required packages

library(data.table)
library(dplyr)
library(genio)

# Read phenotype data
pheno_path <- file.path(results_dir, "simulated_pheno_large.txt")
if (!file.exists(pheno_path)) {
  stop("Phenotype file not found: ", pheno_path)
}
pheno <- fread(pheno_path, colClasses = c("character", "character", "numeric"))

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
    
    # Skip SNPs with insufficient data
    if (n_AA < 5 || n_aa < 5 || n_Aa < 5) {
      results[j,] <- list(snp, NA, NA, NA, NA)
      next
    }
    
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
cat("Processing chromosome", chr, "...\n")
cat("PLINK base:", plink_base, "\n")

if (!all(file.exists(paste0(plink_base, c(".bed", ".bim", ".fam"))))) {
  stop("PLINK files not found for: ", plink_base)
}

plink_data <- read_plink(plink_base)
geno_df <- as.data.frame(t(plink_data$X))
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