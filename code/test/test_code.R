library(data.table)
library(dplyr)
library(BEDMatrix)
## ---- Sample Selection ----

fam_file <- "/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024/ukb22418_c1_b0_v2_s487945.fam"
geno_ids <- fread(fam_file, 
                  header = FALSE, 
                  col.names = c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype"))
geno_ids <- geno_ids[Sex %in% c(1, 2), .(IID)]


pheno_file <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv"
pheno <- fread(pheno_file)

# Filter phenotype data
pheno_filtered <- pheno %>%
  select(
    eid,                           # Participant ID
    age = "21022-0.0",             # Age at recruitment
    sex = "31-0.0",                # Sex (1=Male, 0=Female)
    pc1 = "22009-0.1",             # Genetic PC1
    pc2 = "22009-0.2",             # Genetic PC2
    pc3 = "22009-0.3",             # Genetic PC3
    pc4 = "22009-0.4",             # Genetic PC4
    pc5 = "22009-0.5",             # Genetic PC5
    pc6 = "22009-0.6",             # Genetic PC6
    pc7 = "22009-0.7",             # Genetic PC7
    pc8 = "22009-0.8",             # Genetic PC8
    pc9 = "22009-0.9",             # Genetic PC9
    pc10 = "22009-0.10",           # Genetic PC10
    used_in_pca = "22020-0.0",     # 1 = Used in PCA
    outlier = "22027-0.0",         # 0 = Non-outlier
    ethnic = "22006-0.0"           # Ethnic background
  ) %>%
  filter(
    ethnic == 1,                   # White British
    used_in_pca == 1,            # Used in PCA analysis
    eid %in% geno_ids$IID
  )

# Random sampling
set.seed(1729)
num_samples <- 200000
sample_indices <- sample(nrow(pheno_filtered), num_samples)
sample_ids <- pheno_filtered[sample_indices, .(eid)]

# Create ID file
ID = as.data.frame(sample_ids$eid)
colnames(ID) = "FID"
ID$"IID" = ID$FID
ID_file <- file.path("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom", "sample_ID.txt" )
write.table(ID,  ID_file,
            quote = FALSE, 
            col.names = TRUE,
            row.names = FALSE, 
            sep = "\t")

# Verification: Count unique IDs
cat("Unique IDs in final file:", nrow(ID), "\n")
## ---- PLINK2 filter (per-CHR) ----
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2" 
bed_bim_dir <- "/data2/morgante_lab/data/ukbiobank/genotypes/array"
pheno_dir   <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom"
fam_dir <- "/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024"

maf <- 0.01
mac <- 5
geno <- 0.1
hwe <- 10E-5

threads_plink <- 20

for (CHR in 1:22) {
  bed <- file.path(bed_bim_dir, sprintf("ukb22418_c%d_b0_v2_s488243.bed", CHR))
  bim <- file.path(bed_bim_dir, sprintf("ukb22418_c%d_b0_v2_s488243.bim", CHR))
  fam <- file.path(fam_dir,     sprintf("ukb22418_c%d_b0_v2_s487945.fam", CHR))
  out_prefix <- file.path(pheno_dir, sprintf("ukb_chr%d_filtered", CHR))
  
  args <- c(
    "--bed", bed,
    "--bim", bim,
    "--fam", fam,
    "--keep", ID_file,
    "--rm-dup", "exclude-all",
    "--snps-only", "just-acgt",
    "--maf", as.character(maf),
    "--mac", as.character(mac),
    "--geno", as.character(geno),
    "--hwe", as.character(hwe),
    "--threads", as.character(threads_plink),
    "--make-bed",
    "--out", out_prefix
  )
  
  cat(sprintf("CHR %d filtering...\n", CHR))
  ret <- system2(plink2_path, args = args)
}
cat("All chromosomes filtered.\n")

## ---- Pheno Simulation ----

# (1) Select random causal SNPs genome-wide --------------------------------
# Combine frequencies across chromosomes
geno_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/ukb_chr"
bim_files <- paste0(geno_pattern, 1:22, "_filtered.bim")
cat("\nBIM files:\n")
print(bim_files)

all_snp <- lapply(bim_files, function(bim_file) {
  fread(bim_file, header = FALSE, col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
}) |>
  data.table::rbindlist()
seed_val <- 1266
num_loci <- 100

# Randomly select SNPs
set.seed(seed_val)
selected_indices <- sample(nrow(all_snp), num_loci)

causal_snps <- all_snp[selected_indices,]
output_causal_snps <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/causal_snps.txt"
# Write causal SNPs to extract
write.table(causal_snps, 
            file = output_causal_snps,
            col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
output_causal_snps_IDs <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/causal_snps_IDs.txt"
write.table(causal_snps$snp, 
            file = output_causal_snps_IDs,
            col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE)
# (2) Simulate effects -----------------------------------------------------
set.seed(seed_val)
a <- rnorm(num_loci, mean = 0, sd = sqrt(2))
num_com_dom <- 100
num_incom_dom <- 0
num_additive <- 0

set.seed(seed_val)
dom_type <- rep(c("complete", "incomplete", "additive"), 
                times = c(num_com_dom, num_incom_dom, num_additive))

d <- numeric(num_loci)
# For complete dominance: random sign * a
complete_idx <- which(dom_type == "complete")
d[complete_idx] <- a[complete_idx] * sample(c(-1, 1), length(complete_idx), replace = TRUE)

# For incomplete dominance: random sign * random factor * a
incomplete_idx <- which(dom_type == "incomplete")
d[incomplete_idx] <- a[incomplete_idx] * 
  sample(c(-1, 1), length(incomplete_idx), replace = TRUE) *
  runif(length(incomplete_idx), 0.2, 0.8)

additive_idx <- which(dom_type == "additive")
d[additive_idx] <- 0

true_effects <- data.table(rsid = causal_snps[[2]], a = a, d = d)
true_effects_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/true_effects.txt"
fwrite(true_effects, true_effects_file, sep = "\t")

# (3) Process chromosomes ---------------------------------------------------
# read causal IDs once (file is single-column rsIDs)
causal_ids <- fread(output_causal_snps_IDs, header = FALSE, col.names = "snp")

total_geno <- NULL

for (CHR in 1:22) {
  cat("Processing chromosome", CHR, "\n")
  
  # 该 chr 的 BIM（原始已过滤的 bfile 的 .bim）
  bim_path0 <- paste0(geno_pattern, CHR, "_filtered.bim")
  if (!file.exists(bim_path0)) {
    cat("  - No BIM for chr", CHR, " → skip\n")
    next
  }
  bim0 <- fread(bim_path0, col.names = c("chr","snp","cm","pos","a1","a2"))
  hit_snps <- intersect(bim0$snp, causal_ids$snp)
  
  if (length(hit_snps) == 0) {
    cat("  - No causal SNPs on chr", CHR, " → skip\n")
    next
  }
  
  # 为该 chr 生成仅含命中位点的 snplist，避免空提取
  snplist_chr <- file.path(pheno_dir, sprintf("causal_chr%d.snplist", CHR))
  fwrite(data.table(SNP = hit_snps), snplist_chr, col.names = FALSE)
  
  out_prefix <- file.path(pheno_dir, paste0("temp_chr", CHR))
  
  # 仅提取该 chr 上命中的因果位点
  system(paste(
    shQuote(plink2_path),
    "--bfile",  shQuote(paste0(geno_pattern, CHR, "_filtered")),
    "--extract", shQuote(snplist_chr),
    "--make-bed",
    "--out",    shQuote(out_prefix)
  ))
  
  fam_file <- paste0(out_prefix, ".fam")
  bim_file <- paste0(out_prefix, ".bim")
  bed_path <- paste0(out_prefix, ".bed")
  
  # 若提取后无文件/无位点，跳过
  if (!file.exists(bed_path) || file.size(bed_path) == 0) {
    cat("  - Extract produced empty BED on chr", CHR, " → skip\n")
    unlink(paste0(out_prefix, c(".bed",".bim",".fam",".log")), force = TRUE)
    next
  }
  
  fam <- fread(fam_file, col.names = c("FID","IID","pat","mat","sex","pheno"))
  if (is.null(total_geno)) {
    total_geno <- data.frame(FID = fam$FID, IID = fam$IID, score = 0)
  }
  
  cat("  - Loading Genotype Matrix\n")
  geno <- BEDMatrix(bed_path)  # individuals x SNPs
  if (ncol(geno) == 0) {
    cat("  - No SNPs remaining after extract on chr", CHR, " → skip\n")
    unlink(paste0(out_prefix, c(".bed",".bim",".fam",".log")), force = TRUE)
    next
  }
  
  cat("  - Loading BIM file\n")
  bim  <- fread(bim_file, col.names = c("chr","snp","pos_cM","pos","alt","ref"))
  
  # map BIM SNPs to causal list
  match_idx <- match(bim$snp, causal_ids$snp)
  a_vals <- a[match_idx]; a_vals[is.na(a_vals)] <- 0
  d_vals <- d[match_idx]; d_vals[is.na(d_vals)] <- 0
  
  # accumulate score
  for (j in seq_len(ncol(geno))) {
    g <- geno[, j]
    contrib <- dplyr::case_when(
      g == 0 ~ -a_vals[j],
      g == 1 ~  d_vals[j],
      g == 2 ~  a_vals[j],
      TRUE   ~  0
    )
    total_geno$score <- total_geno$score + contrib
  }
  cat("  - Finished chr", CHR, "\n")
  
  # 清理临时文件
  unlink(paste0(out_prefix, c(".bed",".bim",".fam",".log")), force = TRUE)
}

# (4) Final phenotype calculation -------------------------------------------- 
var_g <- var(total_geno$score) 
h2 = 0.5 
var_e <- ((1/h2)-1)*var_g 
set.seed(seed_val) 
total_geno$pheno <- total_geno$score + rnorm(nrow(total_geno), 0, sqrt(var_e))

pheno <- data.frame( FID = as.character(total_geno$FID), 
                     IID = as.character(total_geno$IID), 
                     pheno = total_geno$pheno )

hist(pheno$pheno)

# Save results
output_pheno <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/simulated_pheno.txt"
write.table(pheno,
            output_pheno,
            row.names = FALSE, quote = FALSE, sep = "\t")
output_geno <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/simulated_geno.txt"
write.table(total_geno[, c("FID", "IID", "score")],
            output_geno,
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("Success! Simulated pheno and geno saved \n")


## ---- Split to training and testing sets ----
set.seed(seed_val)
perm <- sample(nrow(ID))
training_ids <- ID[perm[1:180000], ]
testing_ids  <- ID[perm[180001:200000], ]
output_training_ID <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/training_IDs.txt"
fwrite(training_ids, output_training_ID, sep = "\t", col.names = FALSE)
output_testing_ID <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/testing_IDs.txt"
fwrite(testing_ids,  output_testing_ID,  sep = "\t", col.names = FALSE)

setDT(pheno)
training_pheno <- pheno[IID %in% training_ids$IID]
testing_pheno <- pheno[IID%in% testing_ids$IID,]
output_training_pheno <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/training_pheno.txt"
output_testing_pheno <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test/100comdom/pheno/testing_pheno.txt"
fwrite(training_pheno, output_training_pheno, sep = "\t", col.names = TRUE)
fwrite(testing_pheno,  output_testing_pheno,  sep = "\t", col.names = TRUE)
