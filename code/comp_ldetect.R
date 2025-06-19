library(bigsnpr)
library(data.table)

ldetect_chr <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/prs_effects/other/ldetect-data/EUR/fourier_ls-chr22.bed")
n_block0 <- nrow(ldetect_chr)


bim_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/filter/ukb_chr22_filtered_200k.bim"
ukb_snps <- fread(bim_file, select = c(1, 2, 4), 
                  col.names = c("chr", "snp", "pos"))
setkey(ukb_snps, pos)
ldetect_chr = as.data.table(ldetect_chr)
ldetect_chr[, block_id := .I]

ukb_snps[, ldetect_block := NA_integer_]

for (i in 1:nrow(ldetect_chr)) {
  ukb_snps[
    pos >= ldetect_chr$start[i] & 
      pos < ldetect_chr$stop[i],
    ldetect_block := i
  ]
}

table(ukb_snps$ldetect_block, useNA = "always")

