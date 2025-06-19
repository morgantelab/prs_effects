#------------Test For One Chromosome Only-------------------------

## Load Packgaes
library(data.table)
library(dplyr)

results_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
fam_dir <- "/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024"
pheno_file <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv"

geno_ids <- fread(paste0(fam_dir, "/ukb22418_c1_b0_v2_s487945.fam"), 
                  header = FALSE, 
                  col.names = c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype"))
geno_ids <- geno_ids[Sex %in% c(1, 2), .(IID)]

# Read phenotype data
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
set.seed(123)
sample_indices <- sample(nrow(pheno_filtered), 200000)
sample_ids <- pheno_filtered[sample_indices, .(eid)]

# Create ID file
ID = as.data.frame(sample_ids$eid)
colnames(ID) = "#FID"
ID$"IID" = ID$`#FID`

write.table(ID, file = paste0(results_dir, "/ID_filtered_200k_sim", n, ".txt"), 
            quote = FALSE, 
            col.names = TRUE,
            row.names = FALSE, 
            sep = "\t")

# Verification: Count unique IDs
cat("Unique IDs in final file:", nrow(ID), "\n")
