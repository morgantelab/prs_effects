#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)

score_file <- args[1]
out_prs    <- args[2]
test_ph    <- args[3]
trait      <- args[4]
covar_file <- args[5]
out_corr   <- args[6]

stopifnot(file.exists(score_file))
stopifnot(file.exists(test_ph))
stopifnot(file.exists(covar_file))

# -----------------------------
# Read PRS score
# -----------------------------
dt0 <- fread(score_file)

prs_col <- "SCORE1_SUM"
prs_dt <- dt0[, .(
  IID = as.character(IID),
  PRS = as.numeric(get(prs_col))
)]

fwrite(prs_dt, out_prs, sep = "\t")

# -----------------------------
# Read phenotype and covariates
# -----------------------------
ph_dt <- fread(test_ph)
cv_dt <- fread(covar_file)

ph_dt[, IID := as.character(IID)]
cv_dt[, IID := as.character(IID)]

if (!(trait %in% names(ph_dt))) {
  stop("Trait column not found in phenotype file: ", trait)
}

# keep phenotype
ph_sub <- ph_dt[, .(IID, pheno = as.numeric(get(trait)))]

# identify covariate columns
covars <- setdiff(names(cv_dt), c("FID", "IID"))

# merge phenotype with covariates
ph_work <- merge(ph_sub, cv_dt, by = "IID")
ph_work <- unique(ph_work, by = "IID")

# default: if no covariates, adjusted phenotype = raw phenotype
ph_work[, pheno_adj := pheno]

# -----------------------------
# Phenotype residualization: y ~ C
# -----------------------------
if (length(covars) > 0) {
  model_vars <- c("pheno", covars)
  cc <- complete.cases(ph_work[, ..model_vars])

  if (sum(cc) < 3) {
    stop("Not enough complete observations to fit adjustment model y ~ C")
  }

  rhs <- paste(covars, collapse = " + ")
  fml <- as.formula(paste("pheno ~", rhs))
  fit <- lm(fml, data = ph_work[cc])

  ph_work[cc, pheno_adj := residuals(fit)]
  ph_work[!cc, pheno_adj := NA_real_]
}

ph <- ph_work[, .(IID, pheno_adj)]

# -----------------------------
# Correlation (R and R^2) on adjusted phenotype
# -----------------------------
res <- merge(prs_dt, ph, by = "IID")
cc_corr <- complete.cases(res[, .(PRS, pheno_adj)])

rval <- suppressWarnings(cor(res$pheno_adj[cc_corr], res$PRS[cc_corr]))
r2   <- rval^2

write.table(
  data.frame(corr_beta_r2 = r2, corr_beta_r = rval, N = sum(cc_corr)),
  out_corr,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat(sprintf("[TEST clumped β] N=%d  r=%.6f  R2=%.6f\n", sum(cc_corr), rval, r2))