# summarize_alpha_tables.R
suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
subdirs  <- c("50com50add", "100add", "100comdom")

# helper: read GWAS and try to pull MAF
read_gwas_with_maf <- function(path) {
  if (!file.exists(path)) return(NULL)
  gw <- fread(path)
  # normalize column names we need
  if ("#CHROM" %in% names(gw)) setnames(gw, "#CHROM", "CHR")
  if ("TEST"   %in% names(gw)) gw <- gw[TEST == "ADD"]
  # try possible freq field names -> MAF
  freq_col <- intersect(c("A1_FREQ", "A1FREQ", "A1.FREQ"), names(gw))
  if (length(freq_col) == 1) {
    gw[, MAF := pmin(get(freq_col), 1 - get(freq_col))]
  } else {
    gw[, MAF := as.numeric(NA)]
  }
  # keep minimal set
  keep <- intersect(c("CHR","ID","BETA","MAF"), names(gw))
  gw[, ..keep]
}

for (sd in subdirs) {
  cat("==> Folder:", sd, "\n")
  dir_path <- file.path(base_dir, sd)
  
  # inputs
  true_effects_fp <- file.path(dir_path, "pheno", "true_effects.txt")
  effects_fp      <- file.path(dir_path, "effects_training.txt")
  gwas_fp         <- file.path(dir_path, "gwas_merged.linear")
  
  # -------- read true effects --------
  if (!file.exists(true_effects_fp)) {
    warning("Missing true effects: ", true_effects_fp, " — skipping this folder.")
    next
  }
  te <- fread(true_effects_fp) # expects columns: rsid, a, d
  setnames(te, old = c("rsid","a","d"), new = c("SNP","a_true","d_true"), skip_absent = TRUE)
  
  # -------- read effect estimates (a_hat, d_hat, p_a, p_d) --------
  if (!file.exists(effects_fp)) {
    warning("Missing effects_training.txt in ", dir_path, " — skipping.")
    next
  }
  eff <- fread(effects_fp)
  # unify column names to CHR, SNP, a_hat, d_hat, p_a, p_d
  if ("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")
  # Only keep columns if they exist; fill later as needed
  keep_eff <- intersect(c("CHR","SNP","a_hat","d_hat","p_a","p_d"), names(eff))
  eff <- eff[, ..keep_eff]
  
  # Ensure missing columns exist as NA
  for (nm in c("CHR","SNP","a_hat","d_hat","p_a","p_d")) {
    if (!nm %in% names(eff)) eff[, (nm) := NA]
  }
  
  # -------- read GWAS (beta + MAF if available) --------
  gw <- read_gwas_with_maf(gwas_fp)
  if (is.null(gw)) {
    warning("Missing gwas_merged.linear in ", dir_path, " — will output without beta/MAF.")
    gw <- data.table(CHR=integer(), ID=character(), BETA=numeric(), MAF=numeric())
  }
  
  # -------- merge all --------
  dt <- merge(te, eff, by = "SNP", all.x = TRUE)
  
  # add GWAS beta/MAF on SNP ID
  setnames(gw, "ID", "SNP", skip_absent = TRUE)
  dt <- merge(dt, gw[, .(SNP, CHR_gwas = CHR, beta = BETA, MAF)], by = "SNP", all.x = TRUE)
  
  # Prefer CHR from effects if present; otherwise take from GWAS
  if ("CHR" %in% names(dt)) {
    dt[, CHR_final := CHR]
    dt[is.na(CHR_final) & !is.na(CHR_gwas), CHR_final := CHR_gwas]
    dt[, CHR := CHR_final][, CHR_final := NULL]
  } else {
    setnames(dt, "CHR_gwas", "CHR")
  }
  if ("CHR_gwas" %in% names(dt)) dt[, CHR_gwas := NULL]
  
  # compute alpha (using true a,d and MAF); if MAF missing -> NA
  dt[, alpha := ifelse(is.na(MAF), NA_real_, a_true + d_true * (1 - 2 * MAF))]
  
  # diffs
  if (!("a_hat" %in% names(dt))) dt[, a_hat := NA_real_]
  if (!("d_hat" %in% names(dt))) dt[, d_hat := NA_real_]
  dt[, a_diff := a_true - a_hat]
  dt[, d_diff := d_true - d_hat]
  
  # order and select output columns
  out <- dt[, .(CHR, SNP,
                a_true, d_true, MAF, alpha,
                a_hat, d_hat, a_diff, d_diff,
                p_a, p_d,  # Welch p-values
                beta)]
  
  # write
  out_fp <- file.path(dir_path, "summary_alpha_table.txt")
  fwrite(out, out_fp, sep = "\t")
  cat("   -> wrote:", out_fp, " (", nrow(out), "rows )\n", sep = "")
}

cat("All done.\n")
