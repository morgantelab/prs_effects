library(data.table)

PHENO_CSV  <- "/data2/morgante_lab/data/ukbiobank/download_software/csv/ukb45105.csv"
BCT_sample <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_final_keep.txt"

bct <- fread(BCT_sample, header = TRUE)
setnames(bct, c("FID", "IID"))

get_field_col <- function(cols, field_id) {
  hit <- grep(paste0("^", field_id, "-"), cols, value = TRUE)
  if (length(hit) == 0L) stop("Field ", field_id, " not found in columns.")
  hit[1L]
}

ph_hdr <- fread(PHENO_CSV, nrows = 0)
ph_cols <- names(ph_hdr)

eid_col        <- "eid"
BMI_col        <- get_field_col(ph_cols, "21001")
HEIGHT_col     <- get_field_col(ph_cols, "50")

pheno_needed <- unique(c(
  eid_col, BMI_col, HEIGHT_col
))

pheno <- fread(PHENO_CSV, select = pheno_needed)
setnames(pheno, eid_col, "IID")

# ======== QC BEFORE MERGE WITH BCT ========

setnames(pheno, HEIGHT_col, "height_raw", skip_absent = TRUE)
setnames(pheno, BMI_col,    "BMI_raw",    skip_absent = TRUE)

pheno[, height_raw := as.numeric(height_raw)]
pheno[, BMI_raw    := as.numeric(BMI_raw)]

cat("N total in pheno (subset read):", nrow(pheno), "\n")

cat("\n---- Height raw summary (pheno) ----\n")
print(summary(pheno$height_raw))
cat("Height raw quantiles (0.1%, 1%, 99%, 99.9%):\n")
print(quantile(pheno$height_raw, probs = c(0.001, 0.01, 0.99, 0.999), na.rm = TRUE))

cat("\n---- BMI raw summary (pheno) ----\n")
print(summary(pheno$BMI_raw))
cat("BMI raw quantiles (0.1%, 1%, 99%, 99.9%):\n")
print(quantile(pheno$BMI_raw, probs = c(0.001, 0.01, 0.99, 0.999), na.rm = TRUE))

rankINT <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  qnorm((r - 0.5) / sum(!is.na(r)))
}

# initialize output columns
pheno[, `:=`(height_rint = NA_real_, BMI_rint = NA_real_)]

# RINT independently (no idx needed)
pheno[!is.na(height_raw), height_rint := rankINT(height_raw)]
pheno[!is.na(BMI_raw),    BMI_rint    := rankINT(BMI_raw)]

cat("\nRINT updated independently:\n")
cat("  Height non-missing:", sum(!is.na(pheno$height_raw)), "\n")
cat("  BMI non-missing:",    sum(!is.na(pheno$BMI_raw)), "\n")

# pheno[, mahalanobis := NA_real_]
# pheno[, mahal_outlier := FALSE]
# 
# sub <- pheno[keep_pre_maha == TRUE]
# 
# Y_raw <- as.matrix(sub[, .(height_rint, BMI_rint)])
# ok <- stats::complete.cases(Y_raw)
# Y_raw2 <- Y_raw[ok, , drop = FALSE]
# 
# Y    <- scale(Y_raw2)
# covY <- cov(Y, use = "pairwise.complete.obs")
# 
# D2 <- mahalanobis(
#   x      = Y,
#   center = rep(0, ncol(Y)),
#   cov    = covY
# )
# 
# cutoff <- qchisq(0.99, df = ncol(Y))
# 
# sub_iid_ok <- sub$IID[ok]
# pheno[IID %in% sub_iid_ok, mahalanobis := D2]
# pheno[, mahal_outlier := !is.na(mahalanobis) & mahalanobis > cutoff]
# 
# cat("Step 2 done: Mahalanobis cutoff =", cutoff, "\n")
# cat("Mahalanobis outliers among keep_pre_maha (in pheno):",
#     sum(pheno$mahal_outlier), "\n")

# ======== MERGE QC'ED PHENO WITH BCT ========

samples <- merge(bct, pheno[, .(IID, height_rint, BMI_rint)], by = "IID", all.x = TRUE, sort = FALSE)

# Height phenotype file: keep BCT individuals with height_rint available
height_out <- samples[!is.na(height_rint), .(FID, IID, Height = height_rint)]
cat("\nBCT Height output N:", nrow(height_out), "\n")

fwrite(
  height_out,
  file = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_HEIGHT_RINT.txt",
  sep = "\t",
  na = "NA"
)

# BMI phenotype file: keep BCT individuals with BMI_rint available
bmi_out <- samples[!is.na(BMI_rint), .(FID, IID, BMI = BMI_rint)]
cat("BCT BMI output N:", nrow(bmi_out), "\n")

fwrite(
  bmi_out,
  file = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_RINT.txt",
  sep = "\t",
  na = "NA"
)

# Height keep list (FID IID)
height_keep <- samples[!is.na(height_rint), .(FID, IID)]
fwrite(
  height_keep,
  file = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/keep_HEIGHT.txt",
  sep = "\t",
  col.names = TRUE
)

# BMI keep list (FID IID)
bmi_keep <- samples[!is.na(BMI_rint), .(FID, IID)]
fwrite(
  bmi_keep,
  file = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/keep_BMI.txt",
  sep = "\t",
  col.names = TRUE
)

# Optional: how many in BCT are missing each phenotype
cat("\n---- Missing within BCT after merge ----\n")
cat("BCT total:", nrow(samples), "\n")
cat("BCT missing height:", sum(is.na(samples$height_rint)), "\n")
cat("BCT missing BMI:",    sum(is.na(samples$BMI_rint)), "\n")

# ======== NON-OVERLAP COUNTS AFTER MERGE ========
# 
# n_bct_total <- nrow(bct)
# n_bct_overlap <- sum(!is.na(samples$height_raw) | !is.na(samples$BMI_raw))  # any pheno info present
# n_bct_nonoverlap <- n_bct_total - n_bct_overlap
# 
# n_pheno_total <- nrow(pheno)
# n_pheno_nonoverlap <- sum(!(pheno$IID %in% bct$IID))
# n_pheno_overlap <- n_pheno_total - n_pheno_nonoverlap
# 
# cat("\n---- Overlap summary ----\n")
# cat("BCT total:", n_bct_total, "\n")
# cat("BCT overlap (has any height/BMI):", n_bct_overlap, "\n")
# cat("BCT non-overlap (missing both height & BMI):", n_bct_nonoverlap, "\n")
# cat("Pheno total (subset read):", n_pheno_total, "\n")
# cat("Pheno overlap in BCT:", n_pheno_overlap, "\n")
# cat("Pheno non-overlap (pheno IID not in BCT):", n_pheno_nonoverlap, "\n\n")
# 
# n_total <- nrow(samples)
# 
# n_height_missing <- sum(is.na(samples$height_raw))
# n_BMI_missing    <- sum(is.na(samples$BMI_raw))
# 
# n_both_missing <- sum(
#   is.na(samples$height_raw) & is.na(samples$BMI_raw)
# )
# 
# n_height_only_missing <- sum(
#   is.na(samples$height_raw) & !is.na(samples$BMI_raw)
# )
# 
# n_BMI_only_missing <- sum(
#   !is.na(samples$height_raw) & is.na(samples$BMI_raw)
# )
# 
# n_complete <- sum(
#   !is.na(samples$height_raw) & !is.na(samples$BMI_raw)
# )
# 
# cat("BCT total:", n_total, "\n\n")
# 
# cat("Height missing (any):", n_height_missing, "\n")
# cat("BMI missing (any):", n_BMI_missing, "\n\n")
# 
# cat("Missing BOTH height & BMI:", n_both_missing, "\n")
# cat("Missing ONLY height:", n_height_only_missing, "\n")
# cat("Missing ONLY BMI:", n_BMI_only_missing, "\n\n")
# 
# cat("Complete (height & BMI both present):", n_complete, "\n")
# 
# n_complete <- sum(
#   !is.na(samples$height_raw) & !is.na(samples$BMI_raw)
# )
# 
# # among those, how many are Mahalanobis outliers
# n_outlier_after_complete <- sum(
#   !is.na(samples$height_raw) &
#     !is.na(samples$BMI_raw) &
#     samples$mahal_outlier
# )
# 
# cat("Complete (height & BMI present):", n_complete, "\n")
# cat("Excluded due to Mahalanobis outlier:", n_outlier_after_complete, "\n")
# cat("Remaining after outlier removal:",
#     n_complete - n_outlier_after_complete, "\n")

# ======== OUTPUTS ========

# samples_clean <- samples[
#   keep_pre_maha == TRUE & mahal_outlier == FALSE,
#   .(
#     FID,
#     IID,
#     Height = height_rint,
#     BMI    = BMI_rint
#   )
# ]
# 
# cat("Final clean samples:", nrow(samples_clean), "\n")
# 
# fwrite(
#   samples_clean,
#   file = "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_BMI_HEIGHT.txt",
#   sep = "\t",
#   na = "NA"
# )
