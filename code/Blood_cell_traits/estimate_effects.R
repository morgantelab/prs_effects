#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
})

# -------------------- Args --------------------
args <- commandArgs(trailingOnly = TRUE)
filtered_prefix <- args[1]    # e.g. /.../OUT_DIR/filtered/ukb_chr
snplist_file    <- args[2]    # significant hits
id_list_file    <- args[3]    # FID IID (analysis samples)
pheno_file      <- args[4]    # contains all 16 traits
covar_file      <- args[5]    # covariates used in GWAS
trait_name      <- args[6]    # trait to analyze
output_file     <- args[7]    # final table

# -------------------- Load IDs, phenotype, and covariates --------------------
ids <- fread(id_list_file, header = FALSE, col.names = c("FID","IID"))

ph <- fread(pheno_file)
if (!trait_name %in% names(ph)) {
  stop("Trait not found in phenotype file: ", trait_name)
}
ph <- ph[, .(FID, IID, pheno = get(trait_name))]

covar <- fread(covar_file)
if (!all(c("FID", "IID") %in% names(covar))) {
  stop("Covariate file must contain FID and IID columns: ", covar_file)
}

covar_cols <- setdiff(names(covar), c("FID", "IID"))

ids[,  `:=`(FID = as.character(FID), IID = as.character(IID))]
ph[,   `:=`(FID = as.character(FID),  IID = as.character(IID))]
covar[,`:=`(FID = as.character(FID), IID = as.character(IID))]

sample_df <- merge(ids, ph, by = c("FID", "IID"), all.x = TRUE)
sample_df <- merge(sample_df, covar, by = c("FID", "IID"), all.x = TRUE)

cc <- complete.cases(sample_df[, c("pheno", covar_cols), with = FALSE])
model_df <- sample_df[cc]
fit <- lm(reformulate(covar_cols, response = "pheno"), data = model_df)
model_df[, pheno_resid := residuals(fit)]
pheno_vec <- model_df$pheno_resid
names(pheno_vec) <- model_df$IID

# -------------------- Load SNP list --------------------
snps <- fread(snplist_file, header = FALSE)[[1]]
cat("SNPs to analyze:", length(snps), "\n")

# -------------------- Helper: Welch calculations --------------------
calculate_welch <- function(y, g, snp) {
  AA <- g == 2L
  Aa <- g == 1L
  aa <- g == 0L

  nAA <- sum(AA); nAa <- sum(Aa); naa <- sum(aa)
  if (nAA < 2 || nAa < 2 || naa < 2) return(NULL)

  yAA <- y[AA]; yAa <- y[Aa]; yaa <- y[aa]

  meanAA <- mean(yAA); meanAa <- mean(yAa); meanaa <- mean(yaa)
  varAA  <- var(yAA);  varAa  <- var(yAa);  varaa  <- var(yaa)

  # -------------------- Additive --------------------
  a_hat <- (meanAA - meanaa) / 2
  var_a <- 0.25 * (varAA/nAA + varaa/naa)
  if (!is.finite(var_a) || var_a <= 0) return(NULL)

  df_a <- (varAA/nAA + varaa/naa)^2 /
    ((varAA/nAA)^2/(nAA-1) + (varaa/naa)^2/(naa-1))
  t_a <- a_hat / sqrt(var_a)
  p_a <- 2 * pt(abs(t_a), df = df_a, lower.tail = FALSE)

  # -------------------- Dominance --------------------
  d_hat <- meanAa - (meanAA + meanaa)/2
  var_d <- varAa/nAa + 0.25*(varAA/nAA + varaa/naa)
  if (!is.finite(var_d) || var_d <= 0) return(NULL)

  df_d <- (var_d)^2 / (
    (varAa/nAa)^2/(nAa-1) +
      (varAA/(4*nAA))^2/(nAA-1) +
      (varaa/(4*naa))^2/(naa-1)
  )
  t_d <- d_hat / sqrt(var_d)
  p_d <- 2 * pt(abs(t_d), df = df_d, lower.tail = FALSE)

  data.table(
    SNP = snp,
    a_hat = a_hat,
    d_hat = d_hat,
    p_a   = p_a,
    p_d   = p_d
  )
}

# -------------------- Main Loop: For each chromosome --------------------
final_list <- list()

# -------------------- Impute missing genotypes --------------------
impute_mode <- function(G) {
    for (j in seq_len(ncol(G))) {
      g <- G[, j]
      if (anyNA(g)) {
        tab <- tabulate(g[!is.na(g)] + 1L, nbins = 3)
        g[is.na(g)] <- which.max(tab) - 1L
        G[, j] <- g
      }
    }
    G
}

for (chr in 1:22) {

  cat("\n▶ Chromosome", chr, "\n")

  # filtered bfiles
  base <- sprintf("%s%d_filtered", filtered_prefix, chr)
  bedfile <- paste0(base, ".bed")
  bimfile <- paste0(base, ".bim")
  famfile <- paste0(base, ".fam")

  if (!file.exists(bedfile)) {
    cat("  - missing filtered bed → skip\n")
    next
  }

  # read BIM to find SNPs truly on this chr
  bim <- fread(bimfile, col.names = c("CHR","SNP","cM","BP","ALT","REF"))

  # match SNPs from significant list to this chromosome’s BIM
  idx <- match(snps, bim$SNP)
  keep <- !is.na(idx)

  if (!any(keep)) {
    cat("  - no SNPs on this chromosome\n")
    next
  }

  snp_idx <- idx[keep]
  snp_names <- bim$SNP[snp_idx]

# -------------------- Load genotype matrix --------------------
G <- BEDMatrix(bedfile)[, snp_idx, drop = FALSE]
colnames(G) <- snp_names

 # apply imputation
  G <- impute_mode(G)

  fam <- fread(famfile, col.names = c("FID","IID","pat","mat","sex","pheno"))
  fam_ids <- fam$IID

  # -------------------- Align phenotypes and samples --------------------
  ord <- match(fam_ids, names(pheno_vec))
  y <- pheno_vec[ord]

  ok <- !is.na(y)
  G <- G[ok, , drop = FALSE]
  y <- y[ok]

  cat("  Samples:", length(y), "\n")

  # -------------------- Compute effects --------------------
  out_chr <- lapply(seq_len(ncol(G)), function(j) {
    calculate_welch(y, G[, j], snp_names[j])
  })

  out_chr <- rbindlist(out_chr, use.names = TRUE, fill = TRUE)
  if (nrow(out_chr)) {
    out_chr <- merge(
      out_chr,
      bim[, .(SNP, BP, REF, ALT)],
      by = "SNP",
      all.x = TRUE
    )[, .(
      CHR = chr,
      POS = BP,
      SNP,
      REF,
      ALT,
      a_hat, d_hat, p_a, p_d
    )]

    final_list[[length(final_list) + 1L]] <- out_chr
  }
}

# -------------------- Write Output --------------------
final <- rbindlist(final_list, use.names = TRUE, fill = TRUE)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
fwrite(final, output_file, sep = "\t")

cat("\n✔ Saved results → ", output_file, "\n")
