suppressPackageStartupMessages({
  library(data.table)
})

phi_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_quantnorm.txt"
fam_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/filtered/ukb_chr1_filtered.fam"
out_dir  <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/heritibility/dldsc/phi"

# fam: keep ORIGINAL order
fam <- fread(fam_file, header = FALSE, select = 1:2, col.names = c("FID","IID"))
fam[, rid := .I]
fam[, `:=`(FID = as.character(FID), IID = as.character(IID))]  # <-- 加这一行

tmp <- DT[, .(FID = as.character(get(fid_col)),
              IID = as.character(get(iid_col)),
              PHENO = get(tr))]

DT <- fread(phi_file)

cn <- names(DT)
fid_col <- if ("FID" %in% cn) "FID" else cn[1]
iid_col <- if ("IID" %in% cn) "IID" else cn[2]
traits  <- setdiff(cn, c(fid_col, iid_col))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (tr in traits) {
  
  tmp <- DT[, .(FID = as.character(get(fid_col)),
                IID = as.character(get(iid_col)),
                PHENO = get(tr))]
  
  # IMPORTANT: do NOT reorder fam; join by keys but preserve fam row order via rid
  aligned <- tmp[fam, on = .(FID, IID)]  # keeps fam order (i=fam), but includes rid
  
  # reorder back to original fam line order explicitly
  aligned <- aligned[order(rid), .(FID, IID, PHENO)]
  
  if (nrow(aligned) != nrow(fam)) {
    stop(sprintf("Row count mismatch for trait %s", tr))
  }
  if (any(is.na(aligned$PHENO))) {
    stop(sprintf("Trait %s has NA after matching: phenotype table missing some FID/IID.", tr))
  }
  
  out_file <- file.path(out_dir, paste0(tr, ".phi"))
  fwrite(aligned, file = out_file, sep = "\t", col.names = FALSE)
  cat("Wrote:", out_file, "\n")
}

cat("Done.\n")
