# plot_maf_linear.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
subdirs  <- c("50com50add", "100add", "100comdom")

for (sd in subdirs) {
  message("Processing: ", sd)
  f <- file.path(base_dir, sd, "summary_alpha_table.txt")
  if (!file.exists(f)) {
    warning("Missing file: ", f, " — skipping")
    next
  }
  d <- fread(f)
  
  # compute diffs if not present
  if (!"a_diff" %in% names(d) && all(c("a_true","a_hat") %in% names(d))) {
    d[, a_diff := a_true - a_hat]
  }
  if (!"d_diff" %in% names(d) && all(c("d_true","d_hat") %in% names(d))) {
    d[, d_diff := d_true - d_hat]
  }
  
  # ---- 1) MAF vs a_hat (linear) ----
  p1 <- ggplot(d[!is.na(MAF) & !is.na(a_hat)], aes(MAF, a_hat)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs a_hat"),
         x = "MAF", y = expression(hat(a))) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_a_hat_linear.png"), p1,
         width = 6, height = 4, dpi = 300)
  
  # ---- 2) MAF vs a_diff (linear) ----
  p2 <- ggplot(d[!is.na(MAF) & !is.na(a_diff)], aes(MAF, a_diff)) +
    geom_hline(yintercept = 0, color = "grey60") +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs a_diff"),
         x = "MAF", y = expression(a[true] - hat(a))) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_a_diff_linear.png"), p2,
         width = 6, height = 4, dpi = 300)
  
  # ---- 3) MAF vs d_hat (linear) ----
  p3 <- ggplot(d[!is.na(MAF) & !is.na(d_hat)], aes(MAF, d_hat)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs d_hat"),
         x = "MAF", y = expression(hat(d))) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_d_hat_linear.png"), p3,
         width = 6, height = 4, dpi = 300)
  
  # ---- 4) MAF vs d_diff (linear) ----
  p4 <- ggplot(d[!is.na(MAF) & !is.na(d_diff)], aes(MAF, d_diff)) +
    geom_hline(yintercept = 0, color = "grey60") +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "purple") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs d_diff"),
         x = "MAF", y = expression(d[true] - hat(d))) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_d_diff_linear.png"), p4,
         width = 6, height = 4, dpi = 300)
}


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
subdirs  <- c("50com50add", "100add", "100comdom")

for (sd in subdirs) {
  message("Processing: ", sd)
  f <- file.path(base_dir, sd, "summary_alpha_table.txt")
  if (!file.exists(f)) {
    warning("Missing file: ", f, " — skipping")
    next
  }
  d <- fread(f)
  
  # compute diffs if not present
  if (!"a_diff" %in% names(d) && all(c("a_true","a_hat") %in% names(d))) {
    d[, a_diff := a_true - a_hat]
  }
  if (!"d_diff" %in% names(d) && all(c("d_true","d_hat") %in% names(d))) {
    d[, d_diff := d_true - d_hat]
  }
  
  # ---- 1) MAF vs |a_hat| ----
  p1 <- ggplot(d[!is.na(MAF) & !is.na(a_hat)], aes(MAF, abs(a_hat))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs |a_hat|"),
         x = "MAF", y = expression("|"~hat(a)~"|")) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_abs_a_hat_linear.png"), p1,
         width = 6, height = 4, dpi = 300)
  
  # ---- 2) MAF vs |a_diff| ----
  p2 <- ggplot(d[!is.na(MAF) & !is.na(a_diff)], aes(MAF, abs(a_diff))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs |a_diff|"),
         x = "MAF", y = expression("|"~a[true] - hat(a)~"|")) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_abs_a_diff_linear.png"), p2,
         width = 6, height = 4, dpi = 300)
  
  # ---- 3) MAF vs |d_hat| ----
  p3 <- ggplot(d[!is.na(MAF) & !is.na(d_hat)], aes(MAF, abs(d_hat))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgreen") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs |d_hat|"),
         x = "MAF", y = expression("|"~hat(d)~"|")) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_abs_d_hat_linear.png"), p3,
         width = 6, height = 4, dpi = 300)
  
  # ---- 4) MAF vs |d_diff| ----
  p4 <- ggplot(d[!is.na(MAF) & !is.na(d_diff)], aes(MAF, abs(d_diff))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "purple") +
    scale_x_continuous(limits = c(0, 0.5)) +
    labs(title = paste0(sd, " — MAF vs |d_diff|"),
         x = "MAF", y = expression("|"~d[true] - hat(d)~"|")) +
    theme_bw()
  ggsave(file.path(base_dir, sd, "plot_MAF_vs_abs_d_diff_linear.png"), p4,
         width = 6, height = 4, dpi = 300)
}

