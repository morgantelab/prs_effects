suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

### =========================
### CONFIG
### =========================

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits"
PLOT_DIR <- file.path(OUT_DIR, "Plots_manhattan_pd")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

### =========================
### TRAIT MAP
### =========================

trait_map <- data.table(
  trait_code = c(
    "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
    "Reticulocyte_perc", "HLR_perc", "Platelet_count", "Plateletcrit",
    "PDW", "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
    "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
    "HEIGHT", "BMI"
  ),
  label = c(
    "RBC#", "HGB", "MCV", "RDW", "MSCV",
    "RET%", "HLR%", "PLT#", "PCT", "PDW",
    "WBC#", "LYMPH%", "MONO%", "NEUT%", "EO%", "BASO%",
    "Height", "BMI"
  )
)

### =========================
### FIND FILE
### =========================

get_effect_file <- function(trait) {
  
  f1 <- file.path(
    OUT_DIR, trait, "effects",
    paste0(trait, "_ad.txt")
  )
  
  f2 <- file.path(
    OUT_DIR, "effects",
    paste0(trait, "_ad.txt")
  )
  
  f3 <- file.path(
    OUT_DIR, "HEIGHT_BMI", "effects",
    paste0(trait, "_ad.txt")
  )
  
  if (file.exists(f1)) return(f1)
  if (file.exists(f2)) return(f2)
  if (file.exists(f3)) return(f3)
  
  stop("Missing file for ", trait)
}

### =========================
### LOAD ALL SNP COUNT
### =========================

all_n <- 0L
file_list <- list()

for (i in seq_len(nrow(trait_map))) {
  
  tr <- trait_map$trait_code[i]
  fp <- get_effect_file(tr)
  
  dt <- fread(fp)
  dt[, CHR := as.integer(CHR)]
  dt <- dt[CHR %in% 1:22]
  
  all_n <- all_n + nrow(dt)
  file_list[[tr]] <- fp
}

bonf_thr  <- 0.05 / all_n
bonf_line <- -log10(bonf_thr)

cat("Total SNPs:", all_n, "\n")
cat("Bonferroni threshold:", bonf_thr, "\n")

### =========================
### STYLE
### =========================

my_theme <- theme_classic(base_size = 16) +
  theme(
    text = element_text(family = "sans"),
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(
      size = 19,
      face = "bold",
      hjust = 0.5
    ),
    plot.margin = margin(8, 8, 14, 8)
  )
### =========================
### PREP FUNCTION
### =========================

prep_dt <- function(dt, p_col) {
  
  dt[, CHR := as.integer(CHR)]
  dt[, POS := as.numeric(POS)]
  dt[, p := as.numeric(get(p_col))]
  
  dt <- dt[CHR %in% 1:22]
  dt <- dt[!is.na(p) & p > 0]
  
  dt[p <= 0, p := .Machine$double.xmin]
  
  chr_info <- dt[
    ,
    .(len = max(POS)),
    by = CHR
  ][order(CHR)]
  
  chr_info[
    ,
    start := cumsum(shift(len, fill = 0))
  ]
  
  dt[
    chr_info,
    start := i.start,
    on = "CHR"
  ]
  
  dt[, pos_cum := POS + start]
  
  chr_mid <- dt[
    ,
    .(mid = mean(range(pos_cum))),
    by = CHR
  ][order(CHR)]
  
  dt[, logp := -log10(p)]
  
  list(
    dt = dt,
    chr_mid = chr_mid
  )
}

### =========================
### CREATE INDIVIDUAL PLOTS
### =========================

plot_list <- vector("list", nrow(trait_map))

for (i in seq_len(nrow(trait_map))) {
  
  trait <- trait_map$trait_code[i]
  label <- trait_map$label[i]
  
  dt <- fread(file_list[[trait]])
  
  prep <- prep_dt(dt, "p_d")
  d <- prep$dt
  chr_mid <- prep$chr_mid
  
  d[, grp := CHR %% 2]
  
  p <- ggplot(
    d,
    aes(pos_cum, logp, color = factor(grp))
  ) +
    geom_point(size = 0.8, alpha = 0.85) +
    geom_hline(
      yintercept = bonf_line,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    scale_color_manual(
      values = c(
        "0" = "#2F75B5",  # 深蓝
        "1" = "#D95F59"   # 深红
      ),
      guide = "none"
    )+
    scale_x_continuous(
      breaks = chr_mid$mid,
      labels = chr_mid$CHR,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = label,
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    my_theme
  
  ### Label top significant SNP per chromosome
  
  sig <- d[p < bonf_thr]
  lab_dt <- sig[, .SD[which.min(p)], by = CHR]
  
  if (nrow(lab_dt) > 0) {
    p <- p +
      geom_text_repel(
        data = lab_dt,
        aes(
          x = pos_cum,
          y = logp,
          label = SNP,
          color = factor(grp)
        ),
        size = 4.2,
        family = "sans",
        fontface = "bold",
        max.overlaps = Inf,
        box.padding = 0.45,
        point.padding = 0.25,
        min.segment.length = 0,
        segment.color = "grey35",
        segment.linewidth = 0.35
      )
  }
  
  plot_list[[i]] <- p
  
  cat("✔", trait, "\n")
}

### =========================
### SAVE SIX PLOTS PER PAGE
### TWO PLOTS HORIZONTALLY
### =========================

plots_per_page <- 6
n_pages <- ceiling(length(plot_list) / plots_per_page)

### One combined multi-page PDF
cairo_pdf(
  filename = file.path(PLOT_DIR, "manhattan_pd_all_pages.pdf"),
  width = 20,
  height = 15,
  onefile = TRUE
)
for (pg in seq_len(n_pages)) {
  
  first_plot <- (pg - 1) * plots_per_page + 1
  last_plot <- min(pg * plots_per_page, length(plot_list))
  
  page_plots <- plot_list[first_plot:last_plot]
  
  combined_plot <- wrap_plots(
    page_plots,
    ncol = 2,
    byrow = TRUE
  )
  
  print(combined_plot)
  
}

dev.off()

cat("Saved multi-page PDF:\n")
cat(file.path(PLOT_DIR, "manhattan_pd_all_pages.pdf"), "\n")
cat("Done.\n")
