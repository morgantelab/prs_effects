#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ========= Settings (rep2, p1=0.01, r2=0.5, kb=250) =========
rep             <- 2
plink2_path     <- "/data2/morgante_lab/yuxuans/software/plink2"
pheno_dir       <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/pheno"
results_dir     <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4/test"

pheno_file      <- file.path(pheno_dir, sprintf("simulated_pheno_rep%d.txt", rep))
train_ids_file  <- file.path(pheno_dir, sprintf("training_ID_rep%d.txt", rep))

p1_value        <- 10e-3          # 0.01
r2_value        <- 0.5
clump_kb        <- 250
p_dir_label     <- sprintf("%.2g", p1_value)  # "0.01"

gwas_threads    <- 20
clump_threads   <- 10              # set >0 to pass --memory
keep_intermediates <- FALSE       # set TRUE to keep per-chr files

# ========= Directories & outputs =========
gwas_dir    <- file.path(results_dir, "gwas")
clump_dir   <- file.path(results_dir, "clumping", p_dir_label)
dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)

combined_gwas_path   <- file.path(gwas_dir, sprintf("gwas_rep%d.allchr.glm.linear", rep))
combined_clumps_path <- file.path(clump_dir, sprintf("gwas_clump_rep%d_r2_%.1f.allchr.clumps", rep, r2_value))
combined_snplist_path<- file.path(clump_dir, sprintf("gwas_clump_rep%d_r2_%.1f.allchr.snplist", rep, r2_value))

cat("\n=== GWAS + Clump (combined outputs) ===\n")
cat("rep:", rep, "  p1:", p1_value, "  r2:", r2_value, "  kb:", clump_kb, "\n")
cat("outputs:\n  ", combined_gwas_path, "\n  ", combined_clumps_path, "\n  ", combined_snplist_path, "\n\n")

# ========= Helpers =========
run_plink <- function(args_vec) {
  code <- system2(plink2_path, args = args_vec)
  if (code != 0) warning("PLINK2 exit code ", code)
  invisible(code)
}
add_mem_threads <- function(args_vec, threads, mem_mb = 0) {
  if (threads > 0) args_vec <- c(args_vec, "--threads", as.character(threads))
  if (mem_mb > 0)  args_vec <- c(args_vec, "--memory",  as.character(mem_mb))
  args_vec
}

# ========= Build training phenotype (FID IID pheno) =========
ph <- fread(pheno_file)
stopifnot(all(c("FID","IID","pheno") %in% names(ph)))

tr_ids <- fread(train_ids_file, header = FALSE)
if (ncol(tr_ids) == 1) {
  setnames(tr_ids, "V1", "IID"); tr_ids[, FID := IID]; setcolorder(tr_ids, c("FID","IID"))
} else {
  setnames(tr_ids, c("V1","V2"), c("FID","IID"))
  tr_ids <- tr_ids[, .(FID, IID)]
}
train_pheno <- merge(tr_ids, ph[, .(FID, IID, pheno)], by = c("FID","IID"), all.x = TRUE)
train_pheno_file <- file.path(gwas_dir, sprintf("train_pheno_rep%d.txt", rep))
fwrite(train_pheno, train_pheno_file, sep = "\t")

# ========= 1) GWAS per chr -> combine to one file =========
gwas_files <- character(0)
for (chr in 1:22) {
  bfile_prefix <- file.path(pheno_dir, sprintf("ukb_chr%d_filtered", chr))
  if (!file.exists(paste0(bfile_prefix, ".bim"))) {
    cat("GWAS: skip chr", chr, "(missing BIM)\n"); next
  }
  out_prefix <- file.path(gwas_dir, sprintf("gwas_rep%d_chr%d", rep, chr))
  args <- c(
    "--bfile", bfile_prefix,
    "--pheno", train_pheno_file,
    "--keep",  train_ids_file,
    "--glm", "no-x-sex", "allow-no-covars", "omit-ref",
    "--no-input-missing-phenotype",
    "--pheno-name", "pheno",
    "--out", out_prefix
  )
  args <- add_mem_threads(args, gwas_threads, plink_memory_mb)
  run_plink(args)
  
  gwas_file <- sprintf("%s.pheno.glm.linear", out_prefix)
  if (file.exists(gwas_file)) {
    gwas_files <- c(gwas_files, gwas_file)
  } else {
    warning("Missing GWAS output for chr ", chr)
  }
}

# Combine all GWAS into one table (optional but handy)
if (length(gwas_files)) {
  gwas_list <- lapply(gwas_files, function(p) tryCatch(fread(p), error=function(e) NULL))
  gwas_list <- Filter(Negate(is.null), gwas_list)
  gwas_all  <- rbindlist(gwas_list, fill = TRUE)
  fwrite(gwas_all, combined_gwas_path, sep = "\t")
  cat("Combined GWAS ->", combined_gwas_path, " (rows:", nrow(gwas_all), ")\n")
} else {
  stop("No GWAS outputs found.")
}

# ========= 2) Clump per chr (correct LD) -> combine to one .clumps/.snplist =========
clumps_files <- character(0)
snplist_files <- character(0)

for (chr in 1:22) {
  bfile_prefix <- file.path(pheno_dir, sprintf("ukb_chr%d_filtered", chr))
  gwas_file    <- file.path(gwas_dir, sprintf("gwas_rep%d_chr%d.pheno.glm.linear", rep, chr))
  if (!file.exists(gwas_file)) next
  
  out_prefix <- file.path(clump_dir, sprintf("gwas_clump_rep%d_chr%d_r2_%s",
                                             rep, chr, format(r2_value, nsmall = 1)))
  args <- c(
    "--bfile", bfile_prefix,
    "--clump", gwas_file,
    "--clump-p1", format(p1_value, scientific = FALSE),
    "--clump-r2", format(r2_value, nsmall = 1),
    "--clump-kb", as.character(clump_kb),
    "--keep",    train_ids_file,
    "--out",     out_prefix
  )
  args <- add_mem_threads(args, clump_threads, plink_memory_mb)
  run_plink(args)
  
  clumps_path <- paste0(out_prefix, ".clumps")
  snplist_path <- paste0(out_prefix, ".snplist")
  
  if (file.exists(clumps_path)) clumps_files <- c(clumps_files, clumps_path)
  
  # make per-chr snplist to combine later
  if (file.exists(clumps_path)) {
    dt <- tryCatch(fread(clumps_path), error = function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      id_col <- if ("SNP" %in% names(dt)) {
        "SNP"
      } else if ("ID" %in% names(dt)) {
        "ID"
      } else {
        names(dt)[min(3, ncol(dt))]  # fallback: typical 3rd column
      }
      
      # collect unique IDs as character
      ids <- unique(as.character(dt[[id_col]]))
      ids <- ids[nzchar(ids)]  # drop empty strings
      
      # WRITE as a one-column data.table
      fwrite(data.table(SNP = ids), snplist_path, col.names = FALSE)
      
      snplist_files <- c(snplist_files, snplist_path)
    } else {
      # empty .clumps -> empty snplist
      fwrite(data.table(SNP = character()), snplist_path, col.names = FALSE)
    }
  } else {
    fwrite(data.table(SNP = character()), snplist_path, col.names = FALSE)
  }
}

# Combine all .clumps into one
if (length(clumps_files)) {
  cl_list <- lapply(clumps_files, function(p) tryCatch(fread(p), error=function(e) NULL))
  cl_list <- Filter(Negate(is.null), cl_list)
  cl_all  <- rbindlist(cl_list, fill = TRUE)
  # try to add CHR from files if not present
  if (!("CHR" %in% names(cl_all))) {
    add_chr <- function(p, dt) { dt[, CHR_FILE := as.integer(sub(".*_chr([0-9]+)_.*$", "\\1", basename(p)))] ; dt }
    cl_list2 <- Map(add_chr, clumps_files, cl_list)
    cl_all <- rbindlist(cl_list2, fill = TRUE)
  }
  fwrite(cl_all, combined_clumps_path, sep = "\t")
  cat("Combined clumps ->", combined_clumps_path, " (rows:", nrow(cl_all), ")\n")
} else {
  warning("No .clumps files produced.")
}

# Combine all .snplist into one (unique)
if (length(snplist_files)) {
  vec <- unique(unlist(lapply(snplist_files, function(p) readLines(p, warn = FALSE))))
  vec <- vec[nzchar(vec)]
  fwrite(data.table(SNP = vec), combined_snplist_path, col.names = FALSE)
  cat("Combined snplist ->", combined_snplist_path, " (SNPs:", length(vec), ")\n")
} else {
  warning("No .snplist files produced.")
}

# ========= Cleanup per-chr intermediates if not needed =========
if (!keep_intermediates) {
  # remove per-chr GWAS & clump files, keep the combined ones
  unlink(file.path(gwas_dir, sprintf("gwas_rep%d_chr*", rep)), recursive = FALSE)
  unlink(file.path(clump_dir, sprintf("gwas_clump_rep%d_chr*", rep)), recursive = FALSE)
  cat("Removed per-chr intermediates.\n")
}

cat("\nDone. You now have ONLY the combined outputs:\n  - ", combined_gwas_path,
    "\n  - ", combined_clumps_path, "\n  - ", combined_snplist_path, "\n", sep = "")
