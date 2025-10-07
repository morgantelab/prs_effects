# confusion_matrices_per_folder_VECTOR_STRICT.R
suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
subdirs  <- c("50com50add", "100add", "100comdom")
p_thresh <- 1e-5

safe_div <- function(num, den) ifelse(den > 0, num / den, NA_real_)

for (sd in subdirs) {
  message("\n=== ", sd, " ===")
  dir_path <- file.path(base_dir, sd)
  eff_fp   <- file.path(dir_path, "effects_training.txt")
  te_fp    <- file.path(dir_path, "pheno", "true_effects.txt")
  
  if (!file.exists(eff_fp)) { warning("Missing: ", eff_fp, " — skipping"); next }
  if (!file.exists(te_fp))  { warning("Missing: ", te_fp,  " — skipping");  next }
  
  ## 1) 读取并标准化 ---------------------------------------------------------
  eff_raw <- fread(eff_fp, integer64 = "numeric")
  if ("chromosome" %in% names(eff_raw)) setnames(eff_raw, "chromosome", "CHR")
  
  need_eff <- c("SNP","p_a","p_d")
  if (!all(need_eff %in% names(eff_raw))) {
    stop(sd, ": effects_training.txt must contain SNP, p_a, p_d")
  }
  eff <- eff_raw[, ..need_eff]
  eff[, p_a := as.numeric(p_a)]
  eff[, p_d := as.numeric(p_d)]
  
  # 每个 SNP 只留一行：取最小 p（最显著）
  eff <- eff[, .(
    p_a = suppressWarnings(min(p_a, na.rm = TRUE)),
    p_d = suppressWarnings(min(p_d, na.rm = TRUE))
  ), by = SNP]
  eff[!is.finite(p_a), p_a := NA_real_]
  eff[!is.finite(p_d), p_d := NA_real_]
  
  n_eff_unique <- nrow(eff)
  
  te <- fread(te_fp, integer64 = "numeric")
  setnames(te, old = c("rsid","a","d"), new = c("SNP","a_true","d_true"), skip_absent = TRUE)
  if (!all(c("SNP","a_true","d_true") %in% names(te))) {
    stop(sd, ": true_effects.txt must have rsid/a/d (or SNP/a_true/d_true)")
  }
  
  ## 2) 定义因果集合 ---------------------------------------------------------
  causal_add <- te$SNP                 # 所有因果 = 加性因果
  causal_dom <- te[d_true != 0, SNP]   # 显性因果 = d_true != 0
  
  n_causal_add_total <- length(causal_add)
  n_causal_dom_total <- length(causal_dom)
  
  ## 3) 加性混淆矩阵（用 p_a） -----------------------------------------------
  idx_add_causal      <- match(causal_add, eff$SNP)            # 长度 = #因果
  p_a_causal          <- eff$p_a[idx_add_causal]               # 可能含 NA
  sig_add_causal      <- (!is.na(p_a_causal)) & is.finite(p_a_causal) & (p_a_causal < p_thresh)
  TP_a                <- sum(sig_add_causal, na.rm = TRUE)
  FN_a                <- length(sig_add_causal) - TP_a         # NA 计入 FN
  
  is_noncausal_add    <- !(eff$SNP %chin% causal_add)
  p_a_noncausal       <- eff$p_a[is_noncausal_add]
  sig_add_noncausal   <- (!is.na(p_a_noncausal)) & is.finite(p_a_noncausal) & (p_a_noncausal < p_thresh)
  FP_a                <- sum(sig_add_noncausal, na.rm = TRUE)
  TN_a                <- length(sig_add_noncausal) - FP_a
  
  stopifnot((TP_a + FN_a) == n_causal_add_total)
  stopifnot((FP_a + TN_a) == sum(is_noncausal_add))
  
  cm_add <- data.table(
    Row = c("Causal", "Non-causal"),
    Significant       = c(TP_a, FP_a),
    `Not Significant` = c(FN_a, TN_a)
  )
  
  ## 4) 显性混淆矩阵（用 p_d） -----------------------------------------------
  idx_dom_causal      <- match(causal_dom, eff$SNP)            # 长度 = #显性因果
  p_d_causal          <- eff$p_d[idx_dom_causal]
  sig_dom_causal      <- (!is.na(p_d_causal)) & is.finite(p_d_causal) & (p_d_causal < p_thresh)
  TP_d                <- sum(sig_dom_causal, na.rm = TRUE)
  FN_d                <- length(sig_dom_causal) - TP_d
  
  is_noncausal_dom    <- !(eff$SNP %chin% causal_dom)
  p_d_noncausal       <- eff$p_d[is_noncausal_dom]
  sig_dom_noncausal   <- (!is.na(p_d_noncausal)) & is.finite(p_d_noncausal) & (p_d_noncausal < p_thresh)
  FP_d                <- sum(sig_dom_noncausal, na.rm = TRUE)
  TN_d                <- length(sig_dom_noncausal) - FP_d
  
  stopifnot((TP_d + FN_d) == n_causal_dom_total)
  stopifnot((FP_d + TN_d) == sum(is_noncausal_dom))
  
  cm_dom <- data.table(
    Row = c("Causal", "Non-causal"),
    Significant       = c(TP_d, FP_d),
    `Not Significant` = c(FN_d, TN_d)
  )
  
  ## 5) 指标 ------------------------------------------------------------------
  metrics <- rbindlist(list(
    data.table(
      `Effect type` = "Additive",
      TPR       = safe_div(TP_a, TP_a + FN_a),
      FDR       = safe_div(FP_a, TP_a + FP_a),
      TNR       = safe_div(TN_a, TN_a + FP_a),
      Precision = safe_div(TP_a, TP_a + FP_a),
      Accuracy  = safe_div(TP_a + TN_a, TP_a + TN_a + FP_a + FN_a),
      Power     = safe_div(TP_a, TP_a + FN_a)
    ),
    data.table(
      `Effect type` = "Dominance",
      TPR       = safe_div(TP_d, TP_d + FN_d),
      FDR       = safe_div(FP_d, TP_d + FP_d),
      TNR       = safe_div(TN_d, TN_d + FP_d),
      Precision = safe_div(TP_d, TP_d + FP_d),
      Accuracy  = safe_div(TP_d + TN_d, TP_d + TN_d + FP_d + FN_d),
      Power     = safe_div(TP_d, TP_d + FN_d)
    )
  ), use.names = TRUE)
  
  ## 6) 保存 & 输出 -----------------------------------------------------------
  add_out  <- file.path(dir_path, "confusion_additive.tsv")
  dom_out  <- file.path(dir_path, "confusion_dominance.tsv")
  met_out  <- file.path(dir_path, "confusion_metrics.tsv")
  diag_out <- file.path(dir_path, "confusion_diagnostics.txt")
  
  fwrite(cm_add, add_out, sep = "\t")
  fwrite(cm_dom, dom_out, sep = "\t")
  fwrite(metrics, met_out, sep = "\t")
  
  diag_txt <- sprintf(
    paste0(
      "Folder: %s\n",
      "Unique effect SNPs (after collapse) = %d\n",
      "Additive: #causal=%d,  #non-causal-in-effects=%d\n",
      "  -> TP=%d, FN=%d, FP=%d, TN=%d  [checks: TP+FN==causal, FP+TN==noncausal]\n",
      "Dominance: #causal=%d, #non-causal-in-effects=%d\n",
      "  -> TP=%d, FN=%d, FP=%d, TN=%d  [checks: TP+FN==causal, FP+TN==noncausal]\n"
    ),
    sd, n_eff_unique,
    n_causal_add_total, sum(is_noncausal_add),
    TP_a, FN_a, FP_a, TN_a,
    n_causal_dom_total, sum(is_noncausal_dom),
    TP_d, FN_d, FP_d, TN_d
  )
  writeLines(diag_txt, diag_out)
  
  cat("\nAdditive Effects (", sd, ")\n", sep = ""); print(cm_add)
  cat("\nDominant  Effects (", sd, ")\n", sep = ""); print(cm_dom)
  cat("\nMetrics (", sd, ")\n", sep = ""); print(metrics)
  cat("\nDiagnostics written to: ", diag_out, "\n", sep = "")
  cat("\nSaved:\n  ", add_out, "\n  ", dom_out, "\n  ", met_out, "\n", sep = "")
}
