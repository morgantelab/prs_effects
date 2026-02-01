library(data.table)
library(ggplot2)
library(ggrepel)


EFFECT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/effects"

TRAITS <- c(
  "WBC_count", "RBC_count", "Haemoglobin", "MCV", "RDW",
  "Platelet_count", "Plateletcrit", "PDW", "Lymphocyte_perc",
  "Monocyte_perc", "Neutrophill_perc", "Eosinophill_perc",
  "Basophill_perc", "Reticulocyte_perc", "MSCV", "HLR_perc"
)

# NOTE: filename mapping (Basophil vs Basophill)
TRAIT_TO_FILE <- list(
  WBC_count          = "WBC_count_ad.txt",
  RBC_count          = "RBC_count_ad.txt",
  Haemoglobin        = "Haemoglobin_ad.txt",
  MCV                = "MCV_ad.txt",
  RDW                = "RDW_ad.txt",
  Platelet_count     = "Platelet_count_ad.txt",
  Plateletcrit       = "Plateletcrit_ad.txt",
  PDW                = "PDW_ad.txt",
  Lymphocyte_perc    = "Lymphocyte_perc_ad.txt",
  Monocyte_perc      = "Monocyte_perc_ad.txt",
  Neutrophill_perc   = "Neutrophill_perc_ad.txt",
  Eosinophill_perc   = "Eosinophill_perc_ad.txt",
  Basophill_perc     = "Basophill_perc_ad.txt",
  Reticulocyte_perc  = "Reticulocyte_perc_ad.txt",
  MSCV               = "MSCV_ad.txt",
  HLR_perc           = "HLR_perc_ad.txt"
)

manhattan_plot <- function(dt, p_col, trait, out_png) {
  
  data.table::setDT(dt)
  
  # ---- sanitize types (CRITICAL: avoid integer overflow) ----
  dt[, CHR := as.integer(as.character(CHR))]
  dt[, POS := as.numeric(as.character(POS))]
  
  dt <- dt[CHR %in% 1:22]
  dt <- dt[!is.na(CHR) & !is.na(POS)]
  
  # ---- chr lengths + cumulative starts (numeric cumsum) ----
  chr_info <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR][order(CHR)]
  chr_info[, cum_start := cumsum(data.table::shift(as.numeric(chr_len), fill = 0))]
  
  # ---- add cum_start back to dt via join (no merge headaches) ----
  dt[chr_info, cum_start := i.cum_start, on = "CHR"]
  
  # ---- now pos_cum exists safely as numeric ----
  dt[, pos_cum := POS + cum_start]
  
  # ---- use chromosome span midpoints (NOT mean(pos_cum)) ----
  chr_span <- dt[, .(
    chr_start = min(pos_cum, na.rm = TRUE),
    chr_end   = max(pos_cum, na.rm = TRUE)
  ), by = CHR][order(CHR)]
  chr_span[, mid := (chr_start + chr_end) / 2]
  
  # ---- thresholds ----
  n_snps <- nrow(dt)
  bonf_line <- -log10(0.05 / n_snps)
  gw_line   <- -log10(5e-8)
  
  # ---- -log10(p) ----
  dt[, logp := -log10(as.numeric(get(p_col)))]
  dt[, sig_bonf := get(p_col) < (0.05 / n_snps)]
  
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = pos_cum, y = logp)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.75, color = "#4D79A6") +
    ggplot2::geom_hline(yintercept = bonf_line, color = "red",  linetype = "dashed") +
    ggplot2::geom_hline(yintercept = gw_line,   color = "blue", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = chr_span$mid, labels = chr_span$CHR) +
    ggplot2::labs(
      title = paste("Manhattan Plot —", trait, sprintf("(%s)", p_col)),
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title  = ggplot2::element_text(size = 12)
    )
  label_dt <- dt[sig_bonf == TRUE][
    , .SD[which.min(get(p_col))],
    by = CHR
  ]
  p <- p +
    ggrepel::geom_text_repel(
      data = label_dt,
      aes(label = SNP),
      size = 3
    )
  
  ggplot2::ggsave(out_png, p, width = 12, height = 5, dpi = 150)
  cat("✔ Saved:", out_png, "\n")
}


# ======================================================
# MAIN LOOP
# ======================================================
for (trait in TRAITS) {
  
  file_path <- file.path(EFFECT_DIR, TRAIT_TO_FILE[[trait]])
  
  if (!file.exists(file_path)) {
    cat("❌ Missing:", file_path, "\n")
    next
  }
  
  dt <- fread(file_path)
  
  # Check p_a & p_d exist
  if (!("p_a" %in% names(dt)) | !("p_d" %in% names(dt))) {
    stop("File missing p_a or p_d: ", file_path)
  }
  
  plotdir <- file.path(EFFECT_DIR, "manhattan_plots")
  dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
  
  out_pa <- file.path(plotdir, paste0(trait, "_manhattan_p_a.png"))
  out_pd <- file.path(plotdir, paste0(trait, "_manhattan_p_d.png"))
  
  # Generate both plots
  manhattan_plot(dt, "p_a", trait, out_pa)
  manhattan_plot(dt, "p_d", trait, out_pd)
}

cat("\n🎉 Completed Manhattan plots for p_a and p_d.\n")
