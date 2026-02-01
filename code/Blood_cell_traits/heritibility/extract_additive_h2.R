library(data.table)
library(dplyr)
log_dir <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/heritibility/additive"

# list log files (adjust pattern if your logs are named differently)
log_files <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)

if (length(log_files) == 0) {
  stop("No .log files found in: ", log_dir)
}

extract_h2 <- function(f) {
  x <- readLines(f, warn = FALSE)
  
  # find the line like: "Total Observed scale h2: 0.0655 (0.0074)"
  idx <- grep("^Total Observed scale h2:", x)
  if (length(idx) == 0) {
    return(data.frame(
      trait = sub("\\.log$", "", basename(f)),
      h2_obs = NA_real_,
      se = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  line <- x[idx[1]]
  
  # capture two numbers: h2 and SE in parentheses
  m <- regexec("^Total Observed scale h2:\\s*([0-9eE+\\.-]+)\\s*\\(([0-9eE+\\.-]+)\\)", line)
  g <- regmatches(line, m)[[1]]
  
  if (length(g) < 3) {
    return(data.frame(
      trait = sub("\\.log$", "", basename(f)),
      h2_obs = NA_real_,
      se = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    trait = sub("\\.log$", "", basename(f)),
    h2_obs = as.numeric(g[2]),
    se = as.numeric(g[3]),
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, lapply(log_files, extract_h2))
res <- res[order(res$trait), ]

print(res)

# Save table
out_tsv <- file.path(log_dir, "ldsc_additive_h2_observed.tsv")
out_csv <- file.path(log_dir, "ldsc_additive_h2_observed.csv")

write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(res, out_csv, row.names = FALSE)

cat("Saved:\n", out_tsv, "\n", out_csv, "\n")
