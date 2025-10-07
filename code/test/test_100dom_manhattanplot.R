suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/100causal_test"
dir_path <- file.path(base_dir, "100comdom")

eff <- fread(file.path(dir_path, "effects_training.txt"), integer64 = "numeric")
setnames(eff, old = c("chromosome"), new = c("CHR"), skip_absent = TRUE)
te  <- fread(file.path(dir_path, "pheno", "true_effects.txt"), integer64 = "numeric")
setnames(te, old = c("rsid","a","d"), new = c("SNP","a_true","d_true"), skip_absent = TRUE)

dt <- merge(eff[, .(CHR, position, SNP, p_d)], te, by = "SNP", all.x = TRUE)
dt[, position := as.numeric(position)]
dt[, CHR      := as.integer(CHR)]

dt[, is_causal := !is.na(a_true)]
dt[, causal_type := ifelse(is_causal & d_true != 0, "dominant",
                           ifelse(is_causal & d_true == 0, "additive", NA_character_))]

setorder(dt, CHR, position)
chr_sizes <- dt[!is.na(position), .(chr_len = max(position)), by = CHR]
chr_sizes[, chr_len := as.numeric(chr_len)]
chr_sizes[, offset  := shift(cumsum(chr_len), fill = 0)]
dt <- merge(dt, chr_sizes[, .(CHR, offset)], by = "CHR", all.x = TRUE)
dt[, pos_cum := position + as.numeric(offset)]
ticks <- chr_sizes[order(CHR), .(CHR, mid = offset + chr_len/2)]

sig_line <- -log10(10e-5)

p_dom <- ggplot(dt, aes(x = pos_cum, y = -log10(p_d))) +
  geom_point(color = "#B0B0B0", size = 1, alpha = 0.7) +
  geom_point(data = dt[(is_causal == TRUE) & (causal_type == "dominant")],
             aes(color = causal_type), size = 2.5, na.rm = TRUE) +
  geom_hline(yintercept = sig_line, linetype = "dashed", color = "red") +
  scale_color_manual(
    name   = "Causal loci",
    values = c("additive" = "green3", "dominant" = "orange"),
    labels = c("additive" = "Additive", "dominant" = "Dominant"),
    na.translate = FALSE
  ) +
  scale_x_continuous(breaks = ticks$mid, labels = ticks$CHR) +
  labs(title = "100comdom — Dominance Manhattan",
       x = "Chromosome", y = expression(-log[10](p[d]))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line   = element_line(color = "black"),
    panel.grid  = element_blank(),
    legend.position = "right"
  )
p_dom
ggsave(file.path(dir_path, "manhattan_dominance.png"), p_dom, width = 11, height = 4.5, dpi = 300)
