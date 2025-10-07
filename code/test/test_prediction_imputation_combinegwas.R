#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(genio)
})

# =========================
# CONFIG
# =========================
rep <- 2

# Project paths
root        <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4"
gwas_dir    <- file.path(root, "gwas")
pheno_dir   <- file.path(root, "pheno")
out_dir     <- file.path(root, "test")

# Inputs
effects_file   <- file.path(out_dir, "effects_rep2.txt")                    # has a_hat, d_hat, p_a, p_d
gwas_combined  <- file.path(gwas_dir, sprintf("gwas_rep%d.allchr.glm.linear", rep))
geno_pattern   <- file.path(pheno_dir, "ukb_chr")                           # -> ukb_chr{CHR}_filtered
plink_path     <- "/data2/morgante_lab/yuxuans/software/plink2"

# IDs (LD on training set; PRS on test set)
train_ids_file <- file.path(pheno_dir, sprintf("training_ID_rep%d.txt", rep))
test_ids_file  <- file.path(pheno_dir, sprintf("test_ID_rep%d.txt", rep))

# Step 3 thresholds (significance for union)
# Note: 10E-3 == 0.01 (1e-2). If you intend 1e-3, use 1e-3 explicitly.
p_beta_thr <- 10E-3     # GWAS ADD
p_a_thr    <- 10E-3     # scan p_a
p_d_thr    <- 10E-3     # scan p_d

# Step 4 clumping options
clump_by         <- "beta"       # "a" (use p_a) or "beta" (use GWAS P)
clump_p1         <- 10E-3    # index variant p-value threshold in clumping (P1)
clump_r2         <- 0.5
clump_kb         <- 250
threads_plink    <- 10
plink_memory_mb  <- 0         # set > 0 to pass --memory

# Step 5 PRS_scan options
dominance_p_threshold <- 10E-3 # include d_hat when p_d < this

# Outputs
tmp_dir           <- file.path(out_dir, sprintf("tmp_steps3to5_rep%d", rep))
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

union_tbl_path    <- file.path(out_dir, sprintf("signif_union_rep%d.tsv", rep))
union_snplist     <- file.path(out_dir, sprintf("signif_union_rep%d.snplist", rep))
clump_dir         <- file.path(out_dir, sprintf("clumping_step4_%s", clump_by))
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)
clumps_combined   <- file.path(clump_dir, sprintf("union_clump_rep%d_r2_%.1f.allchr.clumps", rep, clump_r2))
clumps_snplist    <- file.path(clump_dir, sprintf("union_clump_rep%d_r2_%.1f.allchr.snplist", rep, clump_r2))

prs_out           <- file.path(out_dir, sprintf("test_prs_scan2_rep%d.txt", rep))

# =========================
# Helpers
# =========================
run_plink <- function(args_vec) {
  if (threads_plink > 0) args_vec <- c(args_vec, "--threads", as.character(threads_plink))
  if (plink_memory_mb > 0) args_vec <- c(args_vec, "--memory",  as.character(plink_memory_mb))
  code <- system2(plink_path, args = args_vec)
  if (code != 0) warning("PLINK2 exit code ", code)
  invisible(code)
}

load_gwas_add <- function() {
  if (file.exists(gwas_combined)) {
    dt <- fread(gwas_combined)
  } else {
    files <- sprintf("%s/gwas_rep%d_chr%d.pheno.glm.linear", gwas_dir, rep, 1:22)
    files <- files[file.exists(files)]
    if (!length(files)) stop("No GWAS files found.")
    lst <- lapply(files, function(p) tryCatch(fread(p), error=function(e) NULL))
    dt <- rbindlist(Filter(Negate(is.null), lst), fill = TRUE)
  }
  if ("#CHROM" %in% names(dt)) setnames(dt, "#CHROM", "CHR")
  if ("ID" %in% names(dt) && !("SNP" %in% names(dt))) setnames(dt, "ID", "SNP")
  req <- c("CHR","SNP","TEST","P")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("GWAS missing columns: ", paste(miss, collapse=", "))
  dt[, CHR := as.integer(CHR)]
  dt[, SNP := as.character(SNP)]
  add <- dt[TEST == "ADD", .(CHR, SNP, P_beta = P)]
  setorder(add, P_beta)
  add <- add[, .SD[1], by = .(CHR, SNP)]
  add
}

# ---- helpers for mode imputation on genotypes (0/1/2) ----
impute_mode_vec <- function(g) {
  if (!anyNA(g)) return(g)
  tab <- tabulate(g[!is.na(g)] + 1L, nbins = 3L)  # counts for genotypes 0,1,2
  mode_val <- which.max(tab) - 1L                 # ties resolve to smallest genotype
  g[is.na(g)] <- mode_val
  g
}

impute_mode_matrix <- function(X) {
  na_cols <- which(colSums(is.na(X)) > 0L)
  if (length(na_cols)) {
    for (j in na_cols) X[, j] <- impute_mode_vec(X[, j])
  }
  X
}

# =========================
# Step 3: UNION of significant hits
# =========================
cat("\n[Step 3] Building union of significant hits...\n")
eff <- fread(effects_file)
if ("chromosome" %in% names(eff)) setnames(eff, "chromosome", "CHR")
if ("#CHROM" %in% names(eff)) setnames(eff, "#CHROM", "CHR")
req_eff <- c("CHR","SNP","p_a","p_d")
miss_eff <- setdiff(req_eff, names(eff))
if (length(miss_eff)) stop("Effects table missing columns: ", paste(miss_eff, collapse=", "))
eff[, `:=`(CHR = as.integer(CHR), SNP = as.character(SNP))]

gwas_add <- load_gwas_add()

sig_beta <- gwas_add[P_beta < p_beta_thr, .(CHR, SNP, P_beta)]
sig_a    <- eff[p_a    < p_a_thr,       .(CHR, SNP, p_a)]
sig_d    <- eff[p_d    < p_d_thr,       .(CHR, SNP, p_d)]

union_dt <- Reduce(function(x, y) merge(x, y, by = c("CHR","SNP"), all = TRUE),
                   list(sig_beta, sig_a, sig_d))
if (nrow(union_dt) == 0) stop("No significant hits in any test, given current thresholds.")

union_dt[, `:=`(pass_beta = !is.na(P_beta),
                pass_a    = !is.na(p_a),
                pass_d    = !is.na(p_d))]
setcolorder(union_dt, c("CHR","SNP","P_beta","p_a","p_d","pass_beta","pass_a","pass_d"))

fwrite(union_dt, union_tbl_path, sep = "\t")
fwrite(union_dt[, .(SNP)], union_snplist, col.names = FALSE)
cat("Union table  -> ", union_tbl_path, "\n", sep = "")
cat("Union snplist-> ", union_snplist,  "\n", sep = "")

# =========================
# Step 4: Clumping on the union set
# =========================
cat("\n[Step 4] Clumping union set (by =", clump_by, ") ...\n", sep = "")

pcol <- if (clump_by == "a") "p_a" else if (clump_by == "beta") "P_beta" else stop("clump_by must be 'a' or 'beta'.")

perchr_reports <- list()
for (chr in 1:22) {
  chr_union <- union_dt[CHR == chr & !is.na(get(pcol)), .(SNP, P = get(pcol))]
  if (!nrow(chr_union)) next
  fp <- file.path(tmp_dir, sprintf("clump_input_chr%d.tsv", chr))
  fwrite(chr_union, fp, sep = "\t")
  perchrs <- list(input = fp, chr = chr)
  perchrs$out_prefix <- file.path(clump_dir, sprintf("union_clump_rep%d_chr%d_r2_%s", rep, chr, format(clump_r2, nsmall = 1)))
  perchrs$geno  <- sprintf("%s%d_filtered", geno_pattern, chr)
  perchrs$keep  <- train_ids_file
  perchrs$kb    <- clump_kb
  perchrs$r2    <- clump_r2
  perchrs$p1    <- clump_p1
  perchrs$pcol  <- "P"
  perchrs$idcol <- "SNP"
  perchr_reports[[length(perchr_reports) + 1]] <- perchrs
}

if (!length(perchr_reports)) stop("No union variants with p-values for chosen clump_by.")

clump_paths <- character(0)
snplist_paths <- character(0)

for (rec in perchr_reports) {
  if (!file.exists(paste0(rec$geno, ".bim"))) next
  args <- c(
    "--bfile", rec$geno,
    "--clump", rec$input,
    "--clump-id-field", rec$idcol,
    "--clump-p-field",  rec$pcol,
    "--clump-p1", format(rec$p1, scientific = FALSE),
    "--clump-r2", format(rec$r2, nsmall = 1),
    "--clump-kb", as.character(rec$kb),
    "--keep",     rec$keep,
    "--out",      rec$out_prefix
  )
  run_plink(args)
  
  clp <- paste0(rec$out_prefix, ".clumps")
  if (file.exists(clp)) {
    clump_paths <- c(clump_paths, clp)
    dtc <- tryCatch(fread(clp), error=function(e) NULL)
    if (!is.null(dtc) && nrow(dtc) > 0) {
      id_col <- if ("SNP" %in% names(dtc)) "SNP" else if ("ID" %in% names(dtc)) "ID" else names(dtc)[min(3, ncol(dtc))]
      snp_ids <- unique(as.character(dtc[[id_col]]))
      sp <- paste0(rec$out_prefix, ".snplist")
      fwrite(data.table(SNP = snp_ids), sp, col.names = FALSE)
      snplist_paths <- c(snplist_paths, sp)
    }
  }
}

if (!length(clump_paths)) stop("No .clumps produced; check inputs/thresholds.")

cl_list <- lapply(clump_paths, function(p) tryCatch(fread(p), error=function(e) NULL))
cl_list <- Filter(Negate(is.null), cl_list)
cl_all  <- rbindlist(cl_list, fill = TRUE)
fwrite(cl_all, clumps_combined, sep = "\t")

all_snps <- unique(unlist(lapply(snplist_paths, function(p) readLines(p, warn = FALSE))))
all_snps <- all_snps[nzchar(all_snps)]
fwrite(data.table(SNP = all_snps), clumps_snplist, col.names = FALSE)

cat("Combined clumps  -> ", clumps_combined, "\n", sep = "")
cat("Combined snplist -> ", clumps_snplist,  "\n", sep = "")

# =========================
# Step 5: PRS_scan (using effects + selected clumped SNPs)
# =========================
cat("\n[Step 5] Computing PRS_scan on test set...\n")

eff2 <- fread(effects_file)
if ("chromosome" %in% names(eff2)) setnames(eff2, "chromosome", "CHR")
if ("#CHROM" %in% names(eff2)) setnames(eff2, "#CHROM", "CHR")
req_eff2 <- c("CHR","SNP","ref_allele","alt_allele","a_hat","d_hat","p_a","p_d")
miss2 <- setdiff(req_eff2, names(eff2))
if (length(miss2)) stop("Effects table missing columns: ", paste(miss2, collapse=", "))
eff2[, `:=`(CHR = as.integer(CHR), SNP = as.character(SNP))]

clumped_ids <- data.table(SNP = all_snps)
clumped_effects <- merge(clumped_ids, eff2, by = "SNP")
if (!nrow(clumped_effects)) stop("No overlap between clumped SNPs and effects.")

test_ids <- fread(test_ids_file, header = FALSE)[[1]]

prs_vec <- numeric(length(test_ids))
names(prs_vec) <- test_ids

for (chr in 1:22) {
  chr_snps <- clumped_effects[CHR == chr,
                              .(SNP, ref_allele, alt_allele, a_hat, d_hat, p_a, p_d)]
  if (!nrow(chr_snps)) next
  
  # Extract test genotypes for these SNPs
  snplist <- file.path(tmp_dir, sprintf("prs_chr%d.snps", chr))
  fwrite(chr_snps[, .(SNP)], snplist, col.names = FALSE)
  
  bfile <- sprintf("%s%d_filtered", geno_pattern, chr)
  out_prefix <- file.path(tmp_dir, sprintf("prs_extract_chr%d", chr))
  args <- c(
    "--bfile", bfile,
    "--keep",  test_ids_file,
    "--extract", snplist,
    "--make-bed",
    "--out", out_prefix
  )
  run_plink(args)
  
  if (!file.exists(paste0(out_prefix, ".bed"))) next
  gd <- genio::read_plink(out_prefix, verbose = FALSE)
  fam_ids <- gd$fam$id
  snp_ids <- gd$bim$id
  X <- t(gd$X)  # individuals x SNPs
  
  # reorder to test_ids
  if (!all(fam_ids == test_ids)) {
    idx <- match(test_ids, fam_ids)
    if (any(is.na(idx))) stop("Mismatch between test_ids and extracted FAM IDs.")
    X <- X[idx, , drop = FALSE]
    fam_ids <- fam_ids[idx]
  }
  
  # Align effect rows to genotype SNP order
  se <- chr_snps[match(snp_ids, chr_snps$SNP)]
  keep <- !is.na(se$SNP)
  if (!any(keep)) next
  X <- X[, keep, drop = FALSE]
  se <- se[keep]
  
  # Impute missing genotypes by SNP mode
  X <- impute_mode_matrix(X)
  
  # Compute PRS contribution (additive +/- dominance by p_d)
  contrib_tot <- numeric(nrow(X))
  for (j in seq_len(ncol(X))) {
    g  <- X[, j]
    a  <- se$a_hat[j]
    d  <- se$d_hat[j]
    pd <- se$p_d[j]
    
    # additive baseline: -a,0,+a for g=0,1,2
    c_add <- numeric(length(g))
    c_add[g == 0L] <- -a
    c_add[g == 1L] <-  0
    c_add[g == 2L] <-  a
    
    if (!is.na(pd) && pd < dominance_p_threshold) {
      c_dom <- numeric(length(g))
      c_dom[g == 1L] <- d
      contrib_tot <- contrib_tot + c_add + c_dom
    } else {
      contrib_tot <- contrib_tot + c_add
    }
  }
  prs_vec <- prs_vec + contrib_tot
}

# Save PRS_scan
fwrite(data.table(IID = test_ids, PRS_scan = prs_vec), prs_out, sep = "\t")
cat("✅ Saved PRS_scan -> ", prs_out, "\n", sep = "")

# Cleanup
unlink(tmp_dir, recursive = TRUE)
cat("\nDone (Steps 3–5).\n")

