library(data.table)

# 你的 log 文件夹
LOG_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/heritibility/dominance"

# 读单个 log 并提取指标
parse_dldsc_log <- function(f) {
  x <- readLines(f, warn = FALSE)
  
  trait <- sub("\\.log$", "", basename(f))
  
  # ---- h2: 先找包含关键字的那一行，再解析两组数字 ----
  h2_line <- x[grep("Total Observed scale h2", x, fixed = TRUE)]
  h2_line <- if (length(h2_line) > 0L) h2_line[1L] else NA_character_
  
  h2 <- h2_se <- NA_real_
  if (!is.na(h2_line)) {
    m <- regexec(
      "Total Observed scale h2\\s*:\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\(\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\)",
      h2_line, perl = TRUE
    )
    hit <- regmatches(h2_line, m)[[1]]
    if (length(hit) >= 3L) {
      h2    <- as.numeric(hit[2])
      h2_se <- as.numeric(hit[3])
    }
  }
  
  # ---- 其他字段：保持你原来的 grab1 写法也行，这里继续用更稳的“找行+解析”也可以 ----
  grab1 <- function(pattern, group = 1) {
    m <- regexec(pattern, x, perl = TRUE)
    hit <- regmatches(x, m)
    hit <- hit[lengths(hit) > 0L]
    if (length(hit) == 0L) return(NA_character_)
    hit[[1]][group + 1L]
  }
  
  intercept    <- as.numeric(grab1("Intercept:\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\(", 1))
  intercept_se <- as.numeric(grab1("Intercept:\\s*[+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?\\s*\\(\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*\\)", 1))
  lambda_gc    <- as.numeric(grab1("Lambda GC:\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", 1))
  mean_chi2    <- as.numeric(grab1("Mean Chi\\^2:\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)", 1))
  ratio        <- grab1("Ratio:\\s*(\\S+)", 1)
  
  data.table(
    trait = trait,
    log_file = f,
    h2 = h2,
    h2_se = h2_se,
    intercept = intercept,
    intercept_se = intercept_se,
    lambda_gc = lambda_gc,
    mean_chi2 = mean_chi2,
    ratio = ratio
  )
}
# 扫描 folder
log_files <- list.files(LOG_DIR, pattern = "\\.log$", full.names = TRUE)
if (length(log_files) == 0L) stop("No .log files found under: ", LOG_DIR)

# 汇总
res <- rbindlist(lapply(log_files, parse_dldsc_log), fill = TRUE)
setorder(res, trait)

print(res)

# 写出
out_tsv <- file.path(LOG_DIR, "dominance_h2_summary.tsv")
fwrite(res, out_tsv, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote:", out_tsv, "\n")
