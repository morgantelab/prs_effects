#!/usr/bin/env Rscript
# combine_prs_scan2.R
# 用法：
#   Rscript combine_prs_scan2.R <OUT_PATH> <PRS_FILE_1> <PRS_FILE_2> ... <PRS_FILE_N>
#
# 说明：
#   每个 PRS_FILE 是形如 .../r2_0.1/test_prs_scan2_rep{rep}.txt 的两列表（IID\tPRS_scan）。
#   本脚本会读取每个文件，给每行加上对应的 r2，并合并输出为
#   OUT_PATH，列名：IID  r2  PRS_scan（制表符分隔）。

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript combine_prs_scan2.R <OUT_PATH> <PRS_FILE_1> ... <PRS_FILE_N>")
}

out_path <- args[1]
in_files <- args[-1]

# 从路径提取 r2：匹配 r2_数字（允许小数）
extract_r2 <- function(p) {
  m <- regexpr("r2_([0-9]+(?:\\.[0-9]+)?)", p, perl = TRUE)
  if (m[1] == -1) return(NA_real_)
  as.numeric(sub("r2_","",regmatches(p, m)[[1]]))
}

parts <- lapply(in_files, function(f){
  if (!file.exists(f)) {
    warning("File not found: ", f)
    return(NULL)
  }
  r2v <- extract_r2(f)
  dt  <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)

  # 兼容列名：若无列名或不同名，尽量取前两列
  if (!all(c("IID","PRS_scan") %in% names(dt))) {
    if (ncol(dt) < 2L) return(NULL)
    setnames(dt, 1:2, c("IID","PRS_scan"))
  }
  if (!("IID" %in% names(dt) && "PRS_scan" %in% names(dt))) return(NULL)

  dt[, r2 := r2v][, .(IID, r2, PRS_scan)]
})

parts <- Filter(Negate(is.null), parts)
if (length(parts) == 0L) {
  stop("No valid input files to combine.")
}

res <- rbindlist(parts, use.names = TRUE, fill = TRUE)
setcolorder(res, c("IID","r2","PRS_scan"))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_path, sep = "\t")
cat("Wrote combined PRS to:", out_path, "with", nrow(res), "rows\n")
