#------------Test For One Chromosome Only-------------------------

## Load Packgaes
library(data.table)
library(dplyr)
## Check the phenotypes
pheno <- fread("/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv")
gc()
### Sex 31-0.0 Age 21022-0.0 First 10 PCAs 22009-0.1 - 22009-0.10 Ethnic 22006-0.0
### Relatedness 22020-0.0 keep only individuals used in PCA
### Outlier due to heterozygosity and/or rate of missing phenotypes 22027-0.0

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
    outlier = "22027-0.0",          # 0 = Non-outlier, overlapped with used_in_pca
    ethnic = "22006-0.0"
  ) %>%
  filter(
    ethnic == 1,
    used_in_pca == 1              # Keep only individuals used in PCA
  ) %>%
  select(-used_in_pca, -outlier)

### 337464 samples left
gc()

ID = as.data.frame(pheno_filtered$eid)
colnames(ID) = "#FID"
ID$"IID" = ID$`#FID`
write.table(ID, "ID_filtered.txt", 
            quote = FALSE, 
            col.names = TRUE,
            row.names = FALSE, 
            sep = "\t")

fam = fread("/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024/ukb22418_c22_b0_v2_s487945.fam")
sum(fam$V1%in%ID$`#FID`)
