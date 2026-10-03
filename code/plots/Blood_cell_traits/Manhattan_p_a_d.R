suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ============================================================
# 1. Paths
# ============================================================
EFFECT_DIR_BLOOD <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/effects"
EFFECT_DIR_HB    <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/HEIGHT_BMI/effects"

OUT_DIR <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/Plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Trait mapping
# ============================================================
trait_map <- data.table(
  trait_code = c(
    "RBC_count", "Haemoglobin", "MCV", "RDW", "MSCV",
    "Reticulocyte_perc", "HLR_perc", "Platelet_count", "Plateletcrit",
    "PDW", "WBC_count", "Lymphocyte_perc", "Monocyte_perc",
    "Neutrophill_perc", "Eosinophill_perc", "Basophill_perc",
    "HEIGHT", "BMI"
  ),
  Phenotype = c(
    "Red Blood Cell Count",
    "Haemoglobin Concentration",
    "Mean Corpuscular Volume",
    "Red Blood Cell Distribution Width",
    "Mean Sphered Cell Volume",
    "Reticulocyte Percentage",
    "High Light Scatter Reticulocyte Percentage",
    "Platelet Count",
    "Plateletcrit",
    "Platelet Distribution Width",
    "White Blood Cell Count",
    "Lymphocyte Percentage",
    "Monocyte Percentage",
    "Neutrophil Percentage",
    "Eosinophil Percentage",
    "Basophil Percentage",
    "Height",
    "Body Mass Index"
  ),
  Trait = c(
    "RBC#",
    "HGB",
    "MCV",
    "RDW",
    "MSCV",
    "RET%",
    "HLR%",
    "PLT#",
    "PCT",
    "PDW",
    "WBC#",
    "LYMPH%",
    "MONO%",
    "NEUT%",
    "EO%",
    "BASO%",
    "HEIGHT",
    "BMI"
  )
)

# ============================================================
# 3. Peak definition
#    Significant SNPs within <= 1 Mb on the same chromosome
#    are merged into one peak
# ============================================================
PEAK_GAP_BP <- 1e6

# ============================================================
# 4. Helper: read and clean one trait file
# ============================================================
read_trait_effects <- function(file_path) {
  if (!file.exists(file_path)) {
    message("Missing file: ", file_path)
    return(NULL)
  }
  
  dt <- fread(file_path)
  
  required_cols <- c("SNP", "CHR", "POS", "p_d")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", basename(file_path), ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dt <- as.data.table(dt)
  
  dt[, SNP := as.character(SNP)]
  dt[, CHR := as.integer(as.character(CHR))]
  dt[, POS := as.numeric(as.character(POS))]
  dt[, p_d := as.numeric(p_d)]
  
  dt <- dt[
    CHR %in% 1:22 &
      !is.na(SNP) &
      !is.na(CHR) &
      !is.na(POS) &
      !is.na(p_d) &
      p_d > 0
  ]
  
  if (nrow(dt) == 0) return(NULL)
  
  setorder(dt, CHR, POS)
  dt
}

# ============================================================
# 5. First pass:
#    count total number of SNP-trait tests across all traits
# ============================================================
TRAIT_TO_FILE <- c(
  WBC_count         = file.path(EFFECT_DIR_BLOOD, "WBC_count_ad.txt"),
  RBC_count         = file.path(EFFECT_DIR_BLOOD, "RBC_count_ad.txt"),
  Haemoglobin       = file.path(EFFECT_DIR_BLOOD, "Haemoglobin_ad.txt"),
  MCV               = file.path(EFFECT_DIR_BLOOD, "MCV_ad.txt"),
  RDW               = file.path(EFFECT_DIR_BLOOD, "RDW_ad.txt"),
  Platelet_count    = file.path(EFFECT_DIR_BLOOD, "Platelet_count_ad.txt"),
  Plateletcrit      = file.path(EFFECT_DIR_BLOOD, "Plateletcrit_ad.txt"),
  PDW               = file.path(EFFECT_DIR_BLOOD, "PDW_ad.txt"),
  Lymphocyte_perc   = file.path(EFFECT_DIR_BLOOD, "Lymphocyte_perc_ad.txt"),
  Monocyte_perc     = file.path(EFFECT_DIR_BLOOD, "Monocyte_perc_ad.txt"),
  Neutrophill_perc  = file.path(EFFECT_DIR_BLOOD, "Neutrophill_perc_ad.txt"),
  Eosinophill_perc  = file.path(EFFECT_DIR_BLOOD, "Eosinophill_perc_ad.txt"),
  Basophill_perc    = file.path(EFFECT_DIR_BLOOD, "Basophill_perc_ad.txt"),
  Reticulocyte_perc = file.path(EFFECT_DIR_BLOOD, "Reticulocyte_perc_ad.txt"),
  MSCV              = file.path(EFFECT_DIR_BLOOD, "MSCV_ad.txt"),
  HLR_perc          = file.path(EFFECT_DIR_BLOOD, "HLR_perc_ad.txt"),
  
  HEIGHT            = file.path(EFFECT_DIR_HB, "HEIGHT_ad.txt"),
  BMI               = file.path(EFFECT_DIR_HB, "BMI_ad.txt")
)

all_trait_data <- list()
total_tests <- 0L

for (trait_code in names(TRAIT_TO_FILE)) {
  file_path <- TRAIT_TO_FILE[[trait_code]]
  dt <- read_trait_effects(file_path)
  
  if (is.null(dt) || nrow(dt) == 0) next
  
  all_trait_data[[trait_code]] <- dt
  total_tests <- total_tests + nrow(dt)
}

if (length(all_trait_data) == 0 || total_tests == 0) {
  stop("No valid dominance scan data found across traits.")
}

global_dom_thresh <- 0.05 / total_tests

cat("Total SNP-trait tests across all traits:", total_tests, "\n")
cat("Global Bonferroni threshold:", format(global_dom_thresh, scientific = TRUE), "\n")

# ============================================================
# 6. Helper: extract significant dominance peaks for one trait
# ============================================================
extract_dominance_peaks <- function(dt, trait_code, dom_thresh, peak_gap_bp = 1e6) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  n_snps_trait <- nrow(dt)
  
  sig <- copy(dt[p_d < dom_thresh])
  if (nrow(sig) == 0) return(NULL)
  
  setorder(sig, CHR, POS, p_d)
  
  # define peaks within chromosome
  sig[, pos_diff := POS - shift(POS), by = CHR]
  sig[, new_peak := fifelse(is.na(pos_diff) | pos_diff > peak_gap_bp, 1L, 0L)]
  sig[, peak_id_chr := cumsum(new_peak), by = CHR]
  sig[, peak_id := paste0("chr", CHR, "_peak", peak_id_chr)]
  
  peak_summary <- sig[, {
    lead_idx <- which.min(p_d)
    
    lead_snp <- SNP[lead_idx]
    lead_pos <- POS[lead_idx]
    lead_p   <- p_d[lead_idx]
    
    peak_start <- min(POS)
    peak_end   <- max(POS)
    peak_width <- peak_end - peak_start
    
    rest_snps <- SNP[-lead_idx]
    rest_snps_txt <- if (length(rest_snps) == 0) "" else paste(rest_snps, collapse = ";")
    
    .(
      CHR = CHR[1],
      Lead_SNP = lead_snp,
      Lead_POS = lead_pos,
      Lead_p_d = lead_p,
      Peak_start = peak_start,
      Peak_end = peak_end,
      Peak_width_bp = peak_width,
      N_sig_SNPs = .N,
      Rest_SNPs = rest_snps_txt
    )
  }, by = peak_id]
  
  peak_summary[, trait_code := trait_code]
  peak_summary[, n_tested_snps_trait := n_snps_trait]
  peak_summary[, global_dominance_threshold := dom_thresh]
  
  peak_summary <- merge(
    peak_summary,
    trait_map[, .(trait_code, Phenotype, Trait)],
    by = "trait_code",
    all.x = TRUE
  )
  
  setcolorder(
    peak_summary,
    c(
      "Trait", "Phenotype", "trait_code",
      "CHR", "Lead_SNP", "Lead_POS", "Lead_p_d",
      "Peak_start", "Peak_end", "Peak_width_bp",
      "N_sig_SNPs", "Rest_SNPs",
      "n_tested_snps_trait", "global_dominance_threshold",
      "peak_id"
    )
  )
  
  setorder(peak_summary, CHR, Lead_POS)
  peak_summary
}

# ============================================================
# 7. Extract significant dominance peaks for all traits
# ============================================================
all_peak_tables <- list()

for (trait_code in names(all_trait_data)) {
  peak_dt <- extract_dominance_peaks(
    dt = all_trait_data[[trait_code]],
    trait_code = trait_code,
    dom_thresh = global_dom_thresh,
    peak_gap_bp = PEAK_GAP_BP
  )
  
  if (!is.null(peak_dt) && nrow(peak_dt) > 0) {
    all_peak_tables[[trait_code]] <- peak_dt
  }
}

if (length(all_peak_tables) == 0) {
  stop("No significant dominance loci were identified under the global Bonferroni threshold.")
}

dominance_loci_table <- rbindlist(all_peak_tables, use.names = TRUE, fill = TRUE)

# order by trait then genomic position
dominance_loci_table[, Trait := factor(
  Trait,
  levels = trait_map$Trait
)]
setorder(dominance_loci_table, Trait, CHR, Lead_POS)

# ============================================================
# 8. Save full table
# ============================================================
out_tsv <- file.path(OUT_DIR, "significant_dominance_loci_table_global_bonferroni.tsv")
out_csv <- file.path(OUT_DIR, "significant_dominance_loci_table_global_bonferroni.csv")

fwrite(dominance_loci_table, out_tsv, sep = "\t", quote = FALSE)
fwrite(dominance_loci_table, out_csv, quote = TRUE)

# ============================================================
# 9. Save compact manuscript table
# ============================================================
manuscript_table <- copy(dominance_loci_table)[
  , .(
    Trait,
    Phenotype,
    CHR,
    Lead_SNP,
    Lead_POS,
    Lead_p_d,
    Peak_start,
    Peak_end,
    Peak_width_bp,
    N_sig_SNPs,
    Rest_SNPs
  )
]

out_manuscript_tsv <- file.path(
  OUT_DIR,
  "significant_dominance_loci_table_for_manuscript_global_bonferroni.tsv"
)
out_manuscript_csv <- file.path(
  OUT_DIR,
  "significant_dominance_loci_table_for_manuscript_global_bonferroni.csv"
)

fwrite(manuscript_table, out_manuscript_tsv, sep = "\t", quote = FALSE)
fwrite(manuscript_table, out_manuscript_csv, quote = TRUE)

# ============================================================
# 10. Save threshold summary
# ============================================================
threshold_summary <- data.table(
  total_snp_trait_tests = total_tests,
  global_bonferroni_threshold = global_dom_thresh,
  peak_gap_bp = PEAK_GAP_BP
)

out_summary <- file.path(OUT_DIR, "dominance_scan_threshold_summary.tsv")
fwrite(threshold_summary, out_summary, sep = "\t", quote = FALSE)

cat("\nSaved files:\n")
cat("  ", out_tsv, "\n")
cat("  ", out_csv, "\n")
cat("  ", out_manuscript_tsv, "\n")
cat("  ", out_manuscript_csv, "\n")
cat("  ", out_summary, "\n")
