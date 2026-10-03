suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

### =========================
### USER CONFIG
### =========================
OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI"

PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

PHENOS <- c("HEIGHT", "BMI")

PHENO_TO_TRAIT <- list(
  HEIGHT = "HEIGHT",
  BMI    = "BMI"
)

PHENO_LABEL_MAP <- c(
  HEIGHT = "Height",
  BMI    = "BMI"
)

### =========================
### HELPER
### =========================
get_effect_dir <- function(out_dir, pheno) {
  d1 <- file.path(out_dir, pheno, "effects")
  d2 <- file.path(out_dir, "effects")
  if (dir.exists(d1)) return(d1)
  if (dir.exists(d2)) return(d2)
  stop(
    "Cannot find effects directory for pheno = ", pheno,
    "\nTried: ", d1,
    "\n   and: ", d2
  )
}

prepare_manhattan_dt <- function(dt, p_col) {
  setDT(dt)
  
  need_cols <- c("CHR", "POS", "SNP", p_col)
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }
  
  dt <- copy(dt)
  
  dt[, CHR := as.integer(as.character(CHR))]
  dt[, POS := as.numeric(as.character(POS))]
  dt[, pval := as.numeric(as.character(get(p_col)))]
  
  dt <- dt[CHR %in% 1:22]
  dt <- dt[!is.na(CHR) & !is.na(POS) & !is.na(pval)]
  
  dt[pval <= 0, pval := .Machine$double.xmin]
  dt[!is.finite(pval), pval := NA_real_]
  dt <- dt[!is.na(pval)]
  
  if (nrow(dt) == 0) return(NULL)
  
  chr_info <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR][order(CHR)]
  chr_info[, cum_start := cumsum(shift(as.numeric(chr_len), fill = 0))]
  
  dt[chr_info, cum_start := i.cum_start, on = "CHR"]
  dt[, pos_cum := POS + cum_start]
  
  chr_span <- dt[, .(
    chr_start = min(pos_cum, na.rm = TRUE),
    chr_end   = max(pos_cum, na.rm = TRUE)
  ), by = CHR][order(CHR)]
  chr_span[, mid := (chr_start + chr_end) / 2]
  
  dt[, logp := -log10(pval)]
  
  list(
    dt = dt,
    chr_span = chr_span
  )
}

### =========================
### GLOBAL BONFERRONI DENOMINATOR
### =========================
all_n_snps <- 0L
file_info <- list()

for (pheno in PHENOS) {
  effect_dir <- get_effect_dir(OUT_DIR, pheno)
  trait <- PHENO_TO_TRAIT[[pheno]]
  file_path <- file.path(effect_dir, paste0(trait, "_ad.txt"))
  
  if (!file.exists(file_path)) {
    warning("Missing: ", file_path)
    next
  }
  
  dt0 <- fread(file_path)
  
  if (!("p_a" %in% names(dt0)) || !("p_d" %in% names(dt0))) {
    stop("File missing p_a or p_d: ", file_path)
  }
  
  dt0[, CHR := as.integer(as.character(CHR))]
  dt0 <- dt0[CHR %in% 1:22]
  n_this <- nrow(dt0)
  
  all_n_snps <- all_n_snps + n_this
  file_info[[pheno]] <- list(
    file_path = file_path,
    n_snps = n_this
  )
}

if (all_n_snps == 0) {
  stop("No valid SNPs found across Height and BMI.")
}

bonf_thr  <- 0.05 / all_n_snps
bonf_line <- -log10(bonf_thr)

cat("Global Bonferroni denominator =", all_n_snps, "\n")
cat("Global Bonferroni threshold   =", bonf_thr, "\n")

### =========================
### COMMON STYLE
### =========================
my_theme <- theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.title = element_blank()
  )

### =========================
### MANHATTAN PLOT FUNCTION
### =========================
manhattan_plot <- function(dt, p_col, out_png, pheno_label) {
  prep <- prepare_manhattan_dt(dt, p_col)
  
  if (is.null(prep) || nrow(prep$dt) == 0) {
    cat("❌ No valid rows for", out_png, "\n")
    return(NULL)
  }
  
  dt_plot <- prep$dt
  chr_span <- prep$chr_span
  
  dt_plot[, sig_bonf := pval < bonf_thr]
  
  # alternating chromosome colors
  dt_plot[, chr_group := ifelse(CHR %% 2 == 0, "even", "odd")]
  
  p <- ggplot(dt_plot, aes(x = pos_cum, y = logp, color = chr_group)) +
    geom_point(size = 0.8, alpha = 0.75) +
    geom_hline(yintercept = bonf_line, linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c("odd" = "#FA7F6F", "even" = "#82B0D2"), guide = "none") +
    scale_x_continuous(
      breaks = chr_span$mid,
      labels = chr_span$CHR,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = pheno_label,
      y = expression(-log[10](italic(P)))
    ) +
    my_theme
  
  label_dt <- dt_plot[sig_bonf == TRUE, .SD[which.min(pval)], by = CHR]
  
  if (nrow(label_dt) > 0) {
    p <- p +
      geom_text_repel(
        data = label_dt,
        aes(label = SNP),
        size = 3.2,
        family = "Arial",
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.size = 0.25
      )
  }
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 12,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  cat("✔ Saved:", out_png, "\n")
}

### =========================
### MAIN LOOP
### =========================
for (pheno in PHENOS) {
  if (is.null(file_info[[pheno]])) next
  
  trait <- PHENO_TO_TRAIT[[pheno]]
  file_path <- file_info[[pheno]]$file_path
  pheno_label <- PHENO_LABEL_MAP[[pheno]]
  
  dt <- fread(file_path)
  
  out_pa <- file.path(PLOT_DIR, paste0(trait, "_manhattan_p_a.png"))
  out_pd <- file.path(PLOT_DIR, paste0(trait, "_manhattan_p_d.png"))
  
  manhattan_plot(dt, "p_a", out_pa, pheno_label)
  manhattan_plot(dt, "p_d", out_pd, pheno_label)
}

cat("Done.\n")
