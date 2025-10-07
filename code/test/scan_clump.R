#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ========= User inputs =========
plink2_path <- "/data2/morgante_lab/yuxuans/software/plink2"
geno_pattern <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno/ukb_chr"
out_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/test"

# effects_file_2 is already in your session; if not, uncomment:
# effects_file_2 <- fread("/data2/.../effects_rep2.txt")

# Clumping parameters
p1 <- 10e-3   # NOTE: 10e-3 == 0.01; use 1e-3 if you intended 0.001
p2 <- 1
r2 <- 0.5
kb <- 250

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ========= 1) Load BIMs (per-chromosome map: CHR/BP -> SNP ID used by PLINK) =========
bim_by_chr <- vector("list", 22)
for (chr in 1:22) {
  bim_file <- sprintf("%s%d_filtered.bim", geno_pattern, chr)
  if (!file.exists(bim_file)) {
    message("Skip chr ", chr, ": missing ", bim_file)
    next
  }
  bim <- fread(bim_file, col.names = c("CHR_raw","SNP_BIM","CM","BP","A1","A2"))
  bim[, CHR := suppressWarnings(as.integer(CHR_raw))]
  bim[, BP  := as.integer(BP)]
  bim_by_chr[[chr]] <- bim[, .(CHR, BP, SNP_BIM)]
}

# ========= 2) Align effects to BIM IDs by CHR/BP (keep alt_allele for --clump-a1-field) =========
effects <- copy(effects_file_2)
req_cols <- c("chromosome","position","p_a","alt_allele")
stopifnot(all(req_cols %in% names(effects)))

setnames(effects, c("chromosome","position"), c("CHR","BP"))
effects[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP))]
effects <- effects[is.finite(p_a) & p_a > 0]

# Join for each chromosome so we only consider SNPs present in that chr's BIM
aligned_all <- rbindlist(lapply(1:22, function(chr) {
  b <- bim_by_chr[[chr]]
  if (is.null(b)) return(NULL)
  # Join on CHR/BP; keep PLINK ID and carry p_a + alt_allele from effects
  tmp <- b[effects[CHR == chr], on = .(CHR, BP)]
  tmp[!is.na(SNP_BIM), .(CHR, BP, SNP = SNP_BIM, p_a, alt_allele)]
}), fill = TRUE)

# Deduplicate SNPs (if any duplicates)—keep smallest p-value
setorder(aligned_all, p_a)
aligned_all <- aligned_all[, .SD[1], by = SNP]
cat("Aligned variants total:", nrow(aligned_all), "\n")

# ========= 3) Write per-chromosome assoc files (SNP, p_a, alt_allele) and clump =========
clumped_paths <- character(0)

for (chr in 1:22) {
  bfile <- sprintf("%s%d_filtered", geno_pattern, chr)
  if (!file.exists(paste0(bfile, ".bim"))) next
  
  assoc_chr <- aligned_all[CHR == chr, .(SNP, p_a, alt_allele)]
  if (nrow(assoc_chr) == 0) {
    message(sprintf("No aligned SNPs for chr %d; skipping", chr))
    next
  }
  
  assoc_path_chr <- file.path(out_dir, sprintf("rep2_for_clump_chr%d.txt", chr))
  fwrite(assoc_chr, assoc_path_chr, sep = "\t")
  message(sprintf("Chr %d: assoc SNPs written: %d -> %s", chr, nrow(assoc_chr), assoc_path_chr))
  
  out_prefix <- file.path(out_dir, sprintf("rep2_clump_chr%d_%dkb_p1_%g", chr, kb, p1))
  cmd <- sprintf(
    '"%s" --bfile "%s" --clump "%s" --clump-snp-field SNP --clump-field p_a --clump-a1-field alt_allele --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out "%s"',
    plink2_path, bfile, assoc_path_chr, p1, p2, r2, kb, out_prefix
  )
  message(cmd)
  status <- system(cmd)
  if (status != 0) warning("PLINK failed on chr ", chr, " (exit code ", status, ")")
  
  fp <- paste0(out_prefix, ".clumped")
  if (file.exists(fp)) clumped_paths <- c(clumped_paths, fp)
}

# ========= 4) Optional: combine all .clumped outputs =========
combine_clumped <- function(paths) {
  if (!length(paths)) return(data.table())
  rbindlist(lapply(paths, function(p) {
    dt <- tryCatch(fread(p), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    dt[, CLUMP_CHR := suppressWarnings(as.integer(sub(".*_chr([0-9]+)_.*$", "\\1", basename(p))))]
    dt
  }), fill = TRUE)
}

clumped_all <- combine_clumped(clumped_paths)
out_combined <- file.path(out_dir, sprintf("rep2_clump_%dkb_p1_%g.allchr.clumped.tsv", kb, p1))
fwrite(clumped_all, out_combined, sep = "\t")
cat("Index SNPs retained (all autosomes):", nrow(clumped_all), "\n")
cat("Combined clumps -> ", out_combined, "\n")


clumped_files <- list.files(out_dir, pattern = "\\.clumps$", full.names = TRUE)
if (length(clumped_files) == 0) {
  stop("No .clumped files found in: ", out_dir)
}

# Read and combine; add CLUMP_CHR parsed from filename like ..._chr12_...
read_one <- function(p) {
  dt <- tryCatch(fread(p), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  chr <- suppressWarnings(as.integer(sub(".*_chr([0-9]+)_.*$", "\\1", basename(p))))
  dt[, CLUMP_CHR := chr]
  dt[, SOURCE_FILE := basename(p)]
  dt
}

combined <- rbindlist(lapply(clumped_files, read_one), fill = TRUE)
if (!nrow(combined)) stop("No readable .clumped files.")

# Optional: sort by chromosome then by P (if present)
if ("P" %in% names(combined)) {
  setorder(combined, CLUMP_CHR, P)
} else if ("P_BEST" %in% names(combined)) {
  # some PLINK outputs may name it differently; adjust if needed
  setorder(combined, CLUMP_CHR, P_BEST)
}

# Write combined file
combined_path <- file.path(out_dir, "combined_clumps_allchr.tsv")
fwrite(combined, combined_path, sep = "\t")
cat("Wrote combined clumps to:", combined_path, "\n")

# Quick summary
cat("Total index SNP rows:", nrow(combined), "\n")
print(combined[, .N, by = CLUMP_CHR][order(CLUMP_CHR)])


