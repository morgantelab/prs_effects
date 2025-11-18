#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

## ---------- Paths ----------

FAM_DIR    <- "/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024"
PHENO_CSV  <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv"
BLOOD_CSV  <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb673826.csv"
ICD10_CSV  <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb674345.csv"

## ---------- Helper: field ID -> column name ----------

get_field_col <- function(cols, field_id) {
  hit <- grep(paste0("^", field_id, "-"), cols, value = TRUE)
  if (length(hit) == 0L) stop("Field ", field_id, " not found in columns.")
  hit[1L]
}

## --------- 0. Genotyped samples from .fam ---------

fam_files <- list.files(FAM_DIR, pattern = "\\.fam$", full.names = TRUE)
if (length(fam_files) == 0L) stop("No .fam files found in FAM_DIR: ", FAM_DIR)

fam_dt <- rbindlist(lapply(fam_files, fread, header = FALSE))[
  , .(FID = V1, IID = V2)
]
fam_dt <- unique(fam_dt)

cat("N fam samples (unique):", nrow(fam_dt), "\n")

## --------- 1. Baseline / QC fields from PHENO_CSV ---------

ph_hdr <- fread(PHENO_CSV, nrows = 0)
ph_cols <- names(ph_hdr)

eid_col        <- "eid"
sex_col        <- get_field_col(ph_cols, "31")     # self-reported sex
gsex_col       <- get_field_col(ph_cols, "22001")  # genetic sex
preg_col       <- get_field_col(ph_cols, "3140")   # pregnant
hetmiss_col    <- get_field_col(ph_cols, "22027")  # outliers for het/missing
ethnic_col     <- get_field_col(ph_cols, "22006")  # genetic ethnic grouping
pca_cohort_col <- get_field_col(ph_cols, "22020")  # used in genetic PCs

age_col        <- get_field_col(ph_cols, "21022")  # age at recruitment
centre_col     <- get_field_col(ph_cols, "54")     # assessment centre
batch_col      <- get_field_col(ph_cols, "22000")  # genotype measurement batch

## top 10 genetic PCs from 22009
pc_all_cols <- grep("^22009-0\\.", ph_cols, value = TRUE)
pc_cols_10  <- sort(pc_all_cols)[1:10]

pheno_needed <- unique(c(
  eid_col, sex_col, gsex_col, preg_col,
  hetmiss_col, ethnic_col, pca_cohort_col,
  age_col, centre_col, batch_col, pc_cols_10
))

pheno <- fread(PHENO_CSV, select = pheno_needed)
setnames(pheno, eid_col, "IID")

## only genotyped samples
pheno <- pheno[IID %in% fam_dt$IID]
cat("N pheno ∩ fam:", nrow(pheno), "\n")

## --------- 2. Blood cell traits: extract 16 baseline traits ---------

pheno_names <- c(
  "WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW", "Platelet_count",
  "Plateletcrit", "PDW", "Lymphocyte_perc", "Monocyte_perc",
  "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
  "Reticulocyte_perc", "MSCV", "HLR_perc"
)

blood_field_ids <- c(
  "30000", # WBC
  "30010", # RBC
  "30020", # Hemoglobin
  "30040", # MCV
  "30070", # RDW
  "30080", # Platelet count
  "30090", # Platelet crit
  "30110", # Platelet distribution width
  "30180", # Lymphocyte percentage
  "30190", # Monocyte percentage
  "30200", # Neutrophil percentage
  "30210", # Eosinophil percentage
  "30220", # Basophil percentage
  "30240", # Reticulocyte percentage
  "30270", # MSCV
  "30290"  # HLR percentage
)
stopifnot(length(pheno_names) == length(blood_field_ids))

blood_trait_cols <- paste0(blood_field_ids, "-0.0")

blood <- fread(BLOOD_CSV, select = c("eid", blood_trait_cols))
setnames(blood, "eid", "IID")

blood <- blood[IID %in% fam_dt$IID]

setnames(blood,
         old = blood_trait_cols,
         new = pheno_names)

## --------- 3. ICD-10 in-patient data: blood-related disease flag ---------

icd_hdr <- fread(ICD10_CSV, nrows = 0)
icd_cols <- names(icd_hdr)

diag_cols <- grep("^41270-", icd_cols, value = TRUE)
if (length(diag_cols) == 0L) {
  stop("No ICD-10 diagnosis columns matching '^41270-' in ICD10_CSV.")
}

icd <- fread(ICD10_CSV, select = c("eid", diag_cols))
setnames(icd, "eid", "IID")

icd_long <- melt(
  icd,
  id.vars      = "IID",
  measure.vars = diag_cols,
  variable.name = "diag_field",
  value.name    = "icd10"
)[!is.na(icd10)]

blood_icd_prefix <- c(
  "C81","C82","C83","C84","C85","C88","C90","C91","C92","C93","C94","C95",
  "D46",
  "D50","D51","D52","D53","D55","D56","D57","D58","D59","D60","D61","D62","D63","D64",
  "B20","B21","B22","B23","B24",
  "N18","Z49","Y841",
  "K70","K71","K72","K73","K74",
  "D45","E831"
)
blood_icd_regex <- paste0("^(", paste(blood_icd_prefix, collapse = "|"), ")")

icd_long[
  ,
  blood_disease := grepl(blood_icd_regex, icd10)
]

bad_eids <- unique(icd_long[blood_disease == TRUE, IID])
icd_flag <- data.table(
  IID = unique(icd$IID),
  blood_disease = unique(icd$IID) %in% bad_eids
)

## --------- 4. Merge everything ---------

samples <- merge(fam_dt, pheno, by = "IID", all.x = FALSE, all.y = FALSE)
samples <- merge(samples, blood[, c("IID", pheno_names), with = FALSE],
                 by = "IID", all.x = TRUE)
samples <- merge(samples, icd_flag, by = "IID", all.x = TRUE)

samples[is.na(blood_disease), blood_disease := FALSE]

samples[
  ,
  `:=`(
    sex             = get(sex_col),
    gsex            = get(gsex_col),
    pregnant_raw    = get(preg_col),
    hetmiss_val     = get(hetmiss_col),
    genetic_ethnic  = get(ethnic_col),
    pca_used        = get(pca_cohort_col),
    age_raw         = get(age_col),
    assessment_ctr  = get(centre_col),
    geno_batch      = get(batch_col)
  )
]

for (k in seq_along(pc_cols_10)) {
  set(samples,
      j = paste0("pc_genetic", k),
      value = samples[[pc_cols_10[k]]])
}

## --------- 5. Define exclusion flags (except missing blood & Mahalanobis) ---------

# White British
samples[
  ,
  white_british := !is.na(genetic_ethnic) & genetic_ethnic == 1
]

# Sex mismatch: only if both present and disagree
samples[
  ,
  sex_mismatch := !is.na(sex) & !is.na(gsex) & sex != gsex
]

# Het/Miss outliers: 22027 != 0
samples[
  ,
  hetmiss_outlier := (!is.na(hetmiss_val) & hetmiss_val != 0)
]

# PCA cohort (22020 == 1)
samples[
  ,
  in_pca_cohort := !is.na(pca_used) & pca_used == 1
]

# Pregnant: only explicit "1" is TRUE; NA (males) = FALSE
samples[
  ,
  pregnant := (!is.na(pregnant_raw) & pregnant_raw == 1)
]

## --------- 6. Missing blood traits (baseline-only) ---------

samples[
  ,
  missing_blood_trait := rowSums(is.na(.SD)) > 0,
  .SDcols = pheno_names
]

## --------- 7. Pre-Mahalanobis composite keep flag ---------

samples[
  ,
  keep_pre_maha := white_british &
    !sex_mismatch &
    !hetmiss_outlier &
    in_pca_cohort &
    !pregnant &
    !blood_disease &
    !missing_blood_trait
]

cat("N keep_pre_maha:", sum(samples$keep_pre_maha), "\n")

## --------- 8. Quantile normalization of 16 traits (within keep_pre_maha) ---------

rankINT <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  qnorm((r - 0.5) / sum(!is.na(r)))
}

norm_dt <- samples[keep_pre_maha == TRUE, c("IID", pheno_names), with = FALSE]

for (nm in pheno_names) {
  norm_dt[[nm]] <- rankINT(norm_dt[[nm]])
}

samples <- merge(
  samples[, !pheno_names, with = FALSE],
  norm_dt,
  by = "IID",
  all.x = TRUE
)

cat("Quantile normalization done for 16 traits (keep_pre_maha subset)\n")

## --------- 9. Mahalanobis distance (only on keep_pre_maha, center=0 after RINT+scale) ---------

sub <- samples[keep_pre_maha == TRUE]

Y_raw <- as.matrix(sub[, ..pheno_names])
Y     <- scale(Y_raw)   # enforce mean ~0, var ~1 for each trait

covY  <- cov(Y)
D2    <- mahalanobis(Y,
                     center = rep(0, ncol(Y)),  # y' Σ^{-1} y with mean 0
                     cov    = covY)

cutoff <- qchisq(0.99, df = ncol(Y))

# store distances + outliers back into samples
samples[, mahalanobis := NA_real_]
samples[keep_pre_maha == TRUE, mahalanobis := D2]

samples[, mahal_outlier := !is.na(mahalanobis) & mahalanobis > cutoff]

cat("Mahalanobis cutoff =", cutoff, "\n")
cat("Mahalanobis outliers among keep_pre_maha:",
    sum(samples$mahal_outlier), "\n")

## --------- 10. Sequential exclusions with counts (your specified order) ---------

N0 <- nrow(samples)
cat("Total samples at start:", N0, "\n")

valid <- rep(TRUE, N0)

step_info <- data.table(
  step_name   = character(),
  excluded    = integer(),
  remaining   = integer()
)

exclude_step <- function(name, condition) {
  idx <- which(valid & condition)
  n_excl <- length(idx)
  valid[idx] <<- FALSE
  n_rem <- sum(valid)
  
  cat(sprintf("Step: %-26s excluded: %6d   remaining: %6d\n",
              name, n_excl, n_rem))
  
  step_info <<- rbind(
    step_info,
    data.table(
      step_name = name,
      excluded  = n_excl,
      remaining = n_rem
    )
  )
}

# Order: White British → Sex mismatch → Het/Miss → PCA → Pregnant → Blood disease → (we also include missing traits) → Mahalanobis

exclude_step("Non-White British",        !samples$white_british)
exclude_step("Sex mismatch",             samples$sex_mismatch)
exclude_step("Het/miss outlier",         samples$hetmiss_outlier)
exclude_step("Not in PCA cohort 22020",  !samples$in_pca_cohort)
exclude_step("Pregnant",                 samples$pregnant)
exclude_step("Blood-related disease",    samples$blood_disease)
exclude_step("Missing blood traits",     samples$missing_blood_trait)
exclude_step("Mahalanobis outlier",      samples$mahal_outlier)

cat("Final remaining (valid) samples:", sum(valid), "\n")

## --------- 11. Final keep flag (use valid) ---------

samples[, keep := valid]
cat("Final N keep:", sum(samples$keep), "\n")

## --------- 12. Save keep list and phenotypes ---------

keep_dt <- samples[keep == TRUE, .(FID, IID)]
out_keep <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_final_keep.txt"
fwrite(keep_dt, out_keep, sep = "\t", col.names = TRUE)

cat("Wrote keep list to:", out_keep, "\n")

pheno_out <- samples[keep == TRUE, c("FID", "IID", pheno_names), with = FALSE]
out_pheno <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"
fwrite(pheno_out, out_pheno, sep = "\t")

cat("Wrote phenotype file to:", out_pheno, "\n")

## --------- 13. Build covariate file for PLINK2 (your style) ---------

dat <- samples[keep == TRUE,
               .(id = IID,
                 sex,
                 age = as.numeric(age_raw),
                 assessment_centre = factor(assessment_ctr),
                 genotype_measurement_batch = factor(geno_batch),
                 pc_genetic1, pc_genetic2, pc_genetic3, pc_genetic4, pc_genetic5,
                 pc_genetic6, pc_genetic7, pc_genetic8, pc_genetic9, pc_genetic10)]

dat[, age2 := age^2]

Z <- model.matrix(
  ~ sex + age + age2 + assessment_centre + genotype_measurement_batch +
    pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 +
    pc_genetic6 + pc_genetic7 + pc_genetic8 + pc_genetic9 + pc_genetic10,
  data = dat
)

Z <- Z[, -1, drop = FALSE]

cols <- which(colnames(Z) %in% c("age", paste0("pc_genetic", 1:10)))
Z[, cols] <- scale(Z[, cols])

Z[, "age2"] <- Z[, "age"]^2

IID <- dat$id
FID <- dat$id
cov.table <- cbind(FID, IID, Z)

out_covar <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates_plink.txt"
write.table(cov.table,
            file      = out_covar,
            quote     = FALSE,
            row.names = FALSE)

cat("Wrote covariate file to:", out_covar, "\n")

## --------- 14. Save step-wise exclusion summary ---------

out_steps <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_QC_step_counts.tsv"
fwrite(step_info, out_steps, sep = "\t")

cat("Wrote QC step summary to:", out_steps, "\n")
