suppressPackageStartupMessages({
  library(data.table)
})

# ====== 用户自定义 ======
base_dir  <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test2"
scenarios <- c("comdom100", "add100", "mix50_50")
n_reps    <- 10
# ========================

all_rows <- list()

for (scen in scenarios) {
  for (rep in seq_len(n_reps)) {
    f_beta <- file.path(base_dir, scen, sprintf("rep%d", rep), "metrics", "corr_beta.txt")
    f_ad   <- file.path(base_dir, scen, sprintf("rep%d", rep), "metrics", "corr_ad.txt")
    
    if (!file.exists(f_beta) || !file.exists(f_ad)) {
      # 缺文件就跳过
      next
    }
    
    # 读入；有些文件可能带标签行（例如 "corr_beta"），统一转数值并过滤
    beta_raw <- tryCatch(fread(f_beta, header = FALSE)[[1]], error = function(e) NA_character_)
    ad_raw   <- tryCatch(fread(f_ad,   header = FALSE)[[1]], error = function(e) NA_character_)
    
    beta_num <- suppressWarnings(as.numeric(beta_raw))
    ad_num   <- suppressWarnings(as.numeric(ad_raw))
    
    # 只保留数值（丢掉标签/空行）
    beta_num <- beta_num[is.finite(beta_num)]
    ad_num   <- ad_num[is.finite(ad_num)]
    
    # 如果这个 rep 中读到了多行数值，取平均；若为空就设 NA
    beta_val <- if (length(beta_num)) mean(beta_num) else NA_real_
    ad_val   <- if (length(ad_num))   mean(ad_num)   else NA_real_
    
    all_rows[[length(all_rows) + 1]] <- data.table(
      scen = scen, rep = rep,
      PRS_beta = beta_val,
      PRS_ad   = ad_val
    )
  }
}

# 合并
res_dt <- if (length(all_rows)) rbindlist(all_rows, use.names = TRUE) else data.table()

# —— 清洗：确保每个 scen×rep 只有一行；若仍有重复，聚合取均值（极少见）
if (nrow(res_dt)) {
  res_byrep <- res_dt[, .(
    PRS_beta = mean(PRS_beta, na.rm = TRUE),
    PRS_ad   = mean(PRS_ad,   na.rm = TRUE)
  ), by = .(scen, rep)]
} else {
  res_byrep <- data.table()
}

scen_map <- c("add100" = "Additive",
              "mix50_50" = "Mixed",
              "comdom100" = "Complete dominant")

plot_dt <- rbind(
  res_byrep[, .(scen = scen_map[scen], rep, metric = "PRS_beta", value = PRS_beta)],
  res_byrep[, .(scen = scen_map[scen], rep, metric = "PRS_a_d", value = PRS_ad)]
)

# 固定 x 轴顺序
plot_dt[, scen := factor(scen, levels = c("Additive", "Mixed", "Complete dominant"))]

# 绘制 boxplot
p <- ggplot(plot_dt, aes(x = scen, y = value, fill = metric)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.size = 0.8) +
  scale_fill_manual(values = c("PRS_beta" = "#377EB8",   # 蓝色
                               "PRS_a_d"  = "#E41A1C")) + # 红色
  labs(x = "Scenario", y = "Correlation",
       fill = "PRS type",
       title = "PRS correlation with phenotype across scenarios") +
  theme_bw(base_size = 14)

print(p)

ggsave(
  filename = file.path(base_dir, "PRS_correlation_boxplot.png"),
  plot = p,
  width = 8, height = 6, dpi = 300
)

# —— 每个 scen 的 mean ± SE
summary_dt <- if (nrow(res_byrep)) {
  res_byrep[, .(
    mean_beta = mean(PRS_beta, na.rm = TRUE),
    se_beta   = sd(PRS_beta,   na.rm = TRUE) / sqrt(sum(is.finite(PRS_beta))),
    mean_ad   = mean(PRS_ad,   na.rm = TRUE),
    se_ad     = sd(PRS_ad,     na.rm = TRUE) / sqrt(sum(is.finite(PRS_ad)))
  ), by = scen][order(scen)]
} else {
  data.table()
}

# 输出到控制台
cat("\n=== 每个 scen × rep 的相关性（干净版） ===\n")
print(res_byrep)

cat("\n=== 每个 scen 的 mean ± SE ===\n")
print(summary_dt)

# 保存结果
out_all  <- file.path(base_dir, "PRS_corr_all_reps.tsv")
out_sum  <- file.path(base_dir, "PRS_corr_summary.tsv")
fwrite(res_byrep, out_all, sep = "\t")
fwrite(summary_dt, out_sum, sep = "\t")
cat("\n✅ 写出：\n", out_all, "\n", out_sum, "\n", sep = "")
