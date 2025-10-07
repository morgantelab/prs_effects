#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(genio) })

## ========= 配置 =========
rep <- 2
root <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4"
pheno_dir  <- file.path(root, "pheno")
gwas_dir   <- file.path(root, "gwas")
out_dir    <- file.path(root, "test")
plink2     <- "/data2/morgante_lab/yuxuans/software/plink2"
threads    <- 10

# 输入
clumps_path     <- file.path(out_dir, "clumping/0.01", sprintf("gwas_clump_rep%d_r2_0.5.allchr.clumps", rep))
gwas_all_linear <- file.path(gwas_dir, sprintf("gwas_rep%d.allchr.glm.linear", rep))  # 若不存在，将自动拼接分染色体文件
test_ids_file   <- file.path(pheno_dir, sprintf("test_ID_rep%d.txt", rep))
geno_prefix     <- file.path(pheno_dir, "ukb_chr")   # => ukb_chr{CHR}_filtered

# 中间文件与输出
tmp_dir         <- file.path(out_dir, sprintf("tmp_prs_rep%d", rep))
score_file      <- file.path(out_dir, sprintf("prs_scorefile_ADD_rep%d.txt", rep))
prs_score_out   <- file.path(out_dir, sprintf("test_prs_score_rep%d.txt", rep))
prs_gwas_out    <- file.path(out_dir, sprintf("test_prs_gwas_rep%d.txt", rep))

dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

## ========= 小工具 =========
load_gwas_add <- function() {
  if (file.exists(gwas_all_linear)) {
    dt <- fread(gwas_all_linear)
  } else {
    files <- sprintf("%s/gwas_rep%d_chr%d.pheno.glm.linear", gwas_dir, rep, 1:22)
    files <- files[file.exists(files)]
    if (!length(files)) stop("找不到 GWAS 输出（无合并文件且无分染色体文件）")
    lst <- lapply(files, fread)
    dt  <- rbindlist(lst, fill = TRUE)
  }
  if ("#CHROM" %in% names(dt)) setnames(dt, "#CHROM", "CHR")
  if ("ID"     %in% names(dt) && !("SNP" %in% names(dt))) setnames(dt, "ID", "SNP")
  req <- c("CHR","SNP","A1","TEST","BETA")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("GWAS 缺少列: ", paste(miss, collapse=", "))
  dt[, CHR := as.integer(CHR)]
  dt[, SNP := as.character(SNP)]
  # 只取 ADD 行；若重复取 P 最小
  add <- dt[TEST == "ADD", .(CHR, SNP, A1, BETA, P = if ("P" %in% names(dt)) P else NA_real_)]
  setorder(add, P)
  add <- add[, .SD[1], by = .(CHR, SNP)]
  add
}

get_bim_alleles <- function(bim) {
  # 兼容不同 genio 版本：alt/ref 或 a1/a2 或 A1/A2
  if (all(c("alt","ref") %in% names(bim))) {
    list(alt = toupper(bim$alt), ref = toupper(bim$ref))
  } else if (all(c("a1","a2") %in% names(bim))) {
    list(alt = toupper(bim$a1), ref = toupper(bim$a2))
  } else if (all(c("A1","A2") %in% names(bim))) {
    list(alt = toupper(bim$A1), ref = toupper(bim$A2))
  } else {
    stop("BIM 里找不到等位基因列（alt/ref 或 a1/a2 或 A1/A2）")
  }
}

comp_base <- c(A="T", T="A", C="G", G="C")
is_comp <- function(a, b) !is.na(comp_base[as.character(a)]) && comp_base[as.character(a)] == b

## ========= 读取 clumps（索引 SNP） =========
stopifnot(file.exists(clumps_path))
cl <- fread(clumps_path)
if ("#CHROM" %in% names(cl)) setnames(cl, "#CHROM", "CHR")
if ("ID" %in% names(cl) && !("SNP" %in% names(cl))) setnames(cl, "ID", "SNP")
cl <- unique(cl[, .(CHR = as.integer(CHR), SNP = as.character(SNP))])

## ========= 读取 GWAS-ADD 并生成 score 文件（SNP A1 BETA） =========
gwas_add <- load_gwas_add()
score_dt <- merge(cl[, .(SNP)], gwas_add[, .(SNP, A1, BETA)], by = "SNP")
if (nrow(score_dt) == 0) stop("clumps 中的 SNP 在 GWAS-ADD 中均未匹配到")
fwrite(score_dt, score_file, sep = "\t", quote = FALSE)
cat("score 文件 SNP 数：", nrow(score_dt), " -> ", score_file, "\n", sep="")

## ========= 运行 PLINK2 --score（每条染色体；合并） =========
sscore_paths <- character(0)
for (chr in 1:22) {
  bfile <- sprintf("%s%d_filtered", geno_prefix, chr)
  if (!file.exists(paste0(bfile, ".bim"))) next
  snps_chr <- cl[CHR == chr, SNP]
  if (!length(snps_chr)) next
  
  snplist <- file.path(tmp_dir, sprintf("chr%d.snplist", chr))
  fwrite(data.table(SNP = snps_chr), snplist, col.names = FALSE)
  
  outp <- file.path(tmp_dir, sprintf("score_rep%d_chr%d", rep, chr))
  args <- c(
    "--bfile", bfile,
    "--keep",  test_ids_file,
    "--extract", snplist,
    "--score", score_file, "1", "2", "3", "header", "cols=+scoresums", "no-mean-imputation",
    "--out", outp,
    "--threads", as.character(threads)
  )
  code <- system2(plink2, args=args)
  if (code != 0) warning("PLINK2 返回码 ", code, "（chr=", chr, "）")
  f <- paste0(outp, ".sscore")
  if (file.exists(f)) sscore_paths <- c(sscore_paths, f)
}
if (!length(sscore_paths)) stop("没有生成任何 .sscore 文件，请检查日志")

# 合并各 chr 的 .sscore 求和为最终 PRS_score
read_sscore <- function(p) {
  dt <- fread(p)
  sc <- grep("^SCORE\\d+_SUM$", names(dt), value = TRUE)
  if (!length(sc)) stop("未找到 SCORE*_SUM 列：", p)
  dt[, .(IID = as.character(IID), SCORE = as.numeric(get(sc[1])))]
}
parts <- lapply(sscore_paths, read_sscore)
all_scores <- rbindlist(parts, use.names=TRUE)
prs_score <- all_scores[, .(PRS_score = sum(SCORE, na.rm = TRUE)), by = IID]
fwrite(prs_score, prs_score_out, sep = "\t")
cat("✅ 已输出 PLINK --score PRS -> ", prs_score_out, "\n", sep="")

## ========= 手工 g×β（与 score 文件完全一致的 SNP/A1/BETA） =========
test_ids <- fread(test_ids_file, header = FALSE)[[1]]
manual <- data.table(IID = as.character(test_ids), PRS_gwas = 0)

for (chr in 1:22) {
  bfile <- sprintf("%s%d_filtered", geno_prefix, chr)
  if (!file.exists(paste0(bfile, ".bim"))) next
  snps_chr <- cl[CHR == chr, SNP]
  if (!length(snps_chr)) next
  
  snplist <- file.path(tmp_dir, sprintf("chr%d.snplist", chr))
  if (!file.exists(snplist)) fwrite(data.table(SNP = snps_chr), snplist, col.names = FALSE)
  
  outp <- file.path(tmp_dir, sprintf("extract_rep%d_chr%d", rep, chr))
  system2(plink2, args = c("--bfile", bfile, "--keep", test_ids_file,
                           "--extract", snplist, "--make-bed", "--out", outp,
                           "--threads", as.character(threads)))
  if (!file.exists(paste0(outp, ".bed"))) next
  
  gp  <- genio::read_plink(outp, verbose = FALSE)
  
  # IDs and matrix
  fam <- as.character(gp$fam$id)
  X   <- t(gp$X)                     # indiv x snp
  X[is.na(X)] <- 0                   # match no-mean-imputation
  
  # Keep test IDs order
  idx <- match(manual$IID, fam)
  if (any(is.na(idx))) stop("FAM 与 test_ID 不一致")
  X <- X[idx, , drop = FALSE]
  
  # BIM as data.table, keep row order aligned to X columns
  bim_dt <- as.data.table(gp$bim)
  # genio::read_plink usually names the SNP ID column 'id'
  if (!("id" %in% names(bim_dt))) stop("BIM 缺少 'id' 列")
  snp <- as.character(bim_dt$id)
  
  # Take the score rows for the SNPs present here, in the same order
  sc_chr <- score_dt[SNP %in% snp, .(SNP, A1, BETA)]
  sc_chr <- sc_chr[match(snp, sc_chr$SNP)]
  keep <- !is.na(sc_chr$SNP)
  if (!any(keep)) next
  
  X     <- X[, keep, drop = FALSE]
  sc_chr<- sc_chr[keep]
  bim_sub <- bim_dt[keep]            # row-subset safely
  
  # Pull alleles from BIM (supports alt/ref or a1/a2 or A1/A2)
  get_bim_alleles <- function(bim) {
    if (all(c("alt","ref") %in% names(bim))) {
      list(alt = toupper(bim$alt), ref = toupper(bim$ref))
    } else if (all(c("a1","a2") %in% names(bim))) {
      list(alt = toupper(bim$a1), ref = toupper(bim$a2))
    } else if (all(c("A1","A2") %in% names(bim))) {
      list(alt = toupper(bim$A1), ref = toupper(bim$A2))
    } else {
      stop("BIM 里找不到等位基因列（alt/ref 或 a1/a2 或 A1/A2）")
    }
  }
  alle <- get_bim_alleles(bim_sub)
  alt <- alle$alt
  ref <- alle$ref
  a1  <- toupper(sc_chr$A1)
  
  # Align to A1 (counted allele). Flip 2-g when A1==REF (or complement).
  comp <- c(A="T", T="A", C="G", G="C")
  is_comp <- function(a,b) !is.na(comp[a]) && comp[a]==b
  
  flip <- logical(length(a1)); ok <- logical(length(a1))
  for (j in seq_along(a1)) {
    if (is.na(a1[j]) || is.na(alt[j]) || is.na(ref[j])) { ok[j] <- FALSE; next }
    if (a1[j]==alt[j]) { flip[j] <- FALSE; ok[j] <- TRUE
    } else if (a1[j]==ref[j]) { flip[j] <- TRUE;  ok[j] <- TRUE
    } else if (is_comp(a1[j], alt[j])) { flip[j] <- FALSE; ok[j] <- TRUE
    } else if (is_comp(a1[j], ref[j])) { flip[j] <- TRUE;  ok[j] <- TRUE
    } else { ok[j] <- FALSE }
  }
  if (!any(ok)) next
  X     <- X[, ok, drop = FALSE]
  flip  <- flip[ok]
  beta  <- sc_chr$BETA[ok]
  
  if (any(flip)) X[, flip] <- 2 - X[, flip, drop = FALSE]
  
  # g × β contribution
  manual[, PRS_gwas := PRS_gwas + as.numeric(X %*% beta)]
}  

fwrite(manual, prs_gwas_out, sep = "\t")
cat("✅ 已输出手工 g×β PRS -> ", prs_gwas_out, "\n", sep="")

## ========= 清理 =========
unlink(tmp_dir, recursive = TRUE)
cat("完成。\n")
