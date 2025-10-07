# gwas_manhattan.R  (robust to non-numeric P)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
subdirs  <- c("50com50add", "100add", "100comdom")

for (sd in subdirs) {
  message("Processing: ", sd)
  dir_path <- file.path(base_dir, sd)
  gwas_fp  <- file.path(dir_path, "gwas_merged.linear")
  te_fp    <- file.path(dir_path, "pheno", "true_effects.txt")
  
  if (!file.exists(gwas_fp)) { warning("Missing: ", gwas_fp, " — skipping"); next }
  if (!file.exists(te_fp))   { warning("Missing: ", te_fp,   " — skipping"); next }
  
  ## --- Read GWAS (avoid integer64; standardize; coerce types) ---
  gw <- fread(gwas_fp, integer64 = "numeric")
  if ("#CHROM" %in% names(gw)) setnames(gw, "#CHROM", "CHR")
  if ("TEST"   %in% names(gw)) gw <- gw[TEST == "ADD"]
  
  req <- c("CHR","POS","ID","P")
  miss <- setdiff(req, names(gw))
  if (length(miss)) stop(sd, ": gwas_merged.linear missing columns: ", paste(miss, collapse=", "))
  
  gw <- gw[, .(CHR, POS, ID, P)]
  # Coerce columns
  gw[, CHR := as.integer(CHR)]
  gw[, POS := as.numeric(POS)]
  # P may be character; coerce safely
  gw[, P_num := suppressWarnings(as.numeric(P))]
  # Replace zeros with a tiny epsilon to avoid -log10(0)
  eps <- 1e-300
  gw[!is.na(P_num) & P_num == 0, P_num := eps]
  # Keep only finite rows
  gw <- gw[!is.na(CHR) & !is.na(POS) & is.finite(P_num)]
  
  if (nrow(gw) == 0) { warning(sd, ": no valid GWAS rows after coercion"); next }
  
  ## --- Read true effects (causal labels) ---
  te <- fread(te_fp, integer64 = "numeric")
  setnames(te, old = c("rsid","a","d"), new = c("SNP","a_true","d_true"), skip_absent = TRUE)
  need_te <- c("SNP","a_true","d_true")
  if (!all(need_te %in% names(te))) stop(sd, ": true_effects.txt must contain rsid/a/d (or SNP/a_true/d_true).")
  
  te[, causal_type := ifelse(d_true == 0, "additive", "dominant")]
  
  ## --- Merge to tag causal SNPs by rsID ---
  setnames(gw, "ID", "SNP")
  dt <- merge(gw, te[, .(SNP, causal_type)], by = "SNP", all.x = TRUE)
  dt[, is_causal := !is.na(causal_type)]
  
  ## --- Cumulative genomic coordinate for Manhattan ---
  setorder(dt, CHR, POS)
  chr_sizes <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR]
  chr_sizes[, chr_len := as.numeric(chr_len)]
  chr_sizes[, offset  := shift(cumsum(chr_len), fill = 0)]
  dt <- merge(dt, chr_sizes[, .(CHR, offset)], by = "CHR", all.x = TRUE)
  dt[, pos_cum := POS + as.numeric(offset)]
  ticks <- chr_sizes[order(CHR), .(CHR, mid = offset + chr_len/2)]
  
  ## --- Plot ---
  sig_line <- -log10(1e-5)
  p <- ggplot(dt, aes(x = pos_cum, y = -log10(P_num))) +
    geom_point(data = dt[is.na(is_causal) | (is_causal == FALSE)],
               color = "#B0B0B0", size = 1, alpha = 0.7) +
    geom_point(data = dt[is_causal == TRUE],
               aes(color = causal_type), size = 2.5, alpha = 0.95, na.rm = TRUE) +
    geom_hline(yintercept = sig_line, linetype = "dashed", color = "red") +
    scale_color_manual(
      name   = "Causal loci",
      values = c("additive" = "green3", "dominant" = "orange"),
      labels = c("additive" = "Additive", "dominant" = "Dominant"),
      na.translate = FALSE
    ) +
    scale_x_continuous(breaks = ticks$mid, labels = ticks$CHR) +
    labs(title = paste0(sd, " — GWAS Manhattan"),
         x = "Chromosome", y = expression(-log[10](p))) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.line   = element_line(color = "black"),
      panel.grid  = element_blank(),
      legend.position = "right"
    )
  
  out_png <- file.path(dir_path, "gwas_manhattan.png")
  ggsave(out_png, p, width = 11, height = 4.5, dpi = 300)
  message("  → saved: ", out_png)
}
