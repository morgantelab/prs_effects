suppressPackageStartupMessages({
  library(data.table)
  library(BEDMatrix)
  library(ggplot2)
})

# --- paths: update these ---
add_only_file   <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/simulated_geno_additive_only.txt"  # FID IID A_additive
eff_with_d_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/true_effects_with_d.tsv"          # rsid a d
filter_dir      <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/filtered"
causal_ids_file <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/causal_snps_IDs.txt"
plink2_path     <- "/data2/morgante_lab/yuxuans/software/plink2"

# --- load precomputed A_i and effects ---
A_dt <- fread(add_only_file)    # columns: FID IID A_additive
setnames(A_dt, c("FID","IID","A"), c("FID","IID","A"), skip_absent=TRUE)
eff  <- fread(eff_with_d_file)  # columns: rsid a d
stopifnot(all(c("rsid","a","d") %in% names(eff)))
setDT(eff); setkey(eff, rsid)

causal_ids <- fread(causal_ids_file, header=FALSE, col.names="snp")

# quintiles from existing A_i (Q1 lowest, Q5 highest)
qnum <- cut(rank(A_dt$A, ties.method="average")/nrow(A_dt),
            breaks=c(0,.2,.4,.6,.8,1), labels=FALSE, include.lowest=TRUE)
mask_Q1 <- (qnum == 1L)
mask_Q5 <- (qnum == 5L)

# quick check: are we actually using d=a?
tol <- 1e-8
same_sign <- sum(sign(eff$d) == sign(eff$a), na.rm=TRUE)
close_val <- sum(abs(eff$d - eff$a) <= tol, na.rm=TRUE)
cat(sprintf("Same sign(d,a): %d / %d (%.2f%%)\n", same_sign, nrow(eff), 100*same_sign/nrow(eff)))
cat(sprintf("Values near-equal (|d-a|<=%g): %d / %d (%.2f%%)\n",
            tol, close_val, nrow(eff), 100*close_val/nrow(eff)))

# If you also want to verify, in Q5 individuals, the heterozygote contribution equals 'a':
bed_paths <- sort(list.files(filter_dir, pattern="^ukb_chr[0-9]+_filtered\\.bed$", full.names=TRUE))
prefixes  <- sub("\\.bed$", "", bed_paths)
stopifnot(length(prefixes) > 0)

obs_list <- list()
for (pref in prefixes) {
  bimf <- paste0(pref,".bim"); bedf <- paste0(pref,".bed"); famf <- paste0(pref,".fam")
  if (!file.exists(bimf) || !file.exists(bedf) || !file.exists(famf)) next
  
  bim  <- fread(bimf, col.names=c("chr","snp","cm","pos","a1","a2"))
  hits <- intersect(bim$snp, causal_ids$snp); if (!length(hits)) next
  
  snplist <- tempfile(fileext=".snplist")
  fwrite(data.table(SNP=hits), snplist, col.names=FALSE)
  outp <- tempfile(pattern="extract_")
  cmd <- sprintf('%s --bfile %s --extract %s --make-bed --out %s',
                 shQuote(plink2_path), shQuote(pref), shQuote(snplist), shQuote(outp))
  ret <- system(cmd)
  if (ret!=0 || !file.exists(paste0(outp,".bed"))) {
    unlink(c(snplist, paste0(outp,c(".bed",".bim",".fam",".log",".nosex"))), force=TRUE); next
  }
  
  fam <- fread(paste0(outp,".fam"), col.names=c("FID","IID","pat","mat","sex","pheno"))
  fam[, `:=`(FID=as.character(FID), IID=as.character(IID))]
  idx <- match(fam$IID, A_dt$IID)
  Q5_chr <- mask_Q5[idx]
  
  G <- BEDMatrix(paste0(outp,".bed"))
  if (is.null(ncol(G)) || ncol(G)==0) {
    unlink(c(snplist, paste0(outp,c(".bed",".bim",".fam",".log",".nosex"))), force=TRUE); next
  }
  
  bim2 <- fread(paste0(outp,".bim"), col.names=c("chr","snp","pos_cM","pos","alt","ref"))
  eff_block <- eff[J(bim2$snp)]
  eff_block[is.na(a), a := 0]; eff_block[is.na(d), d := 0]
  
  # mean contribution among Q5 heterozygotes: if rule is d=a for Q5, this should equal a
  mu_contrib <- rep(NA_real_, ncol(G))
  for (j in seq_len(ncol(G))) {
    g <- G[, j]
    het_Q5 <- which(Q5_chr & g==1L)
    if (!length(het_Q5)) next
    # contribution used in your scoring for hets in Q5:
    # if your scoring uses exactly 'd' in Q5 and you set d=a, the mean below equals a
    mu_contrib[j] <- mean(rep(eff_block$d[j], length(het_Q5)))
  }
  obs_list[[length(obs_list)+1L]] <- data.table(rsid=bim2$snp,
                                                a=eff_block$a, d=eff_block$d,
                                                mu_contrib=mu_contrib)
  
  unlink(c(snplist, paste0(outp,c(".bed",".bim",".fam",".log",".nosex"))), force=TRUE)
}

obs <- rbindlist(obs_list, use.names=TRUE, fill=TRUE)
obs <- obs[!is.na(mu_contrib)]

cat(sprintf("Q5 heterozygote mean contribution vs a:\n  Corr=%.4f, mean(|mu-a|)=%.4g, median(|mu-a|)=%.4g\n",
            suppressWarnings(cor(obs$mu_contrib, obs$a)),
            mean(abs(obs$mu_contrib - obs$a)),
            median(abs(obs$mu_contrib - obs$a))))

# quick plot: expected a vs observed mean contribution in Q5 hets
ggplot(obs, aes(x=a, y=mu_contrib)) +
  geom_point(alpha=.5, size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title="Q5 heterozygotes: observed mean contribution vs a",
       x="a (expected if d=a)", y="observed mean contribution") +
  theme_bw()



suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## If you already have A and S in memory, skip the fread parts below.
A <- fread("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/simulated_geno_additive_only.txt")
S <- fread("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/simulated_geno.txt")
stopifnot(all(A$IID == S$IID))

## Quintiles by additive-only score (Q1 lowest, Q5 highest)
qnum <- cut(rank(A$A_additive, ties.method = "average") / nrow(A),
            breaks = c(0, .2, .4, .6, .8, 1),
            labels = FALSE, include.lowest = TRUE)
qlab <- factor(qnum, levels = 1:5,
               labels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))

## Dominance contribution per individual = final - additive
dom_contrib <- S$score - A$A_additive

## Summarize mean dominance per quintile
dom_summ <- data.table(quintile = qlab, dom = dom_contrib)[
  , .(mean_dom = mean(dom, na.rm = TRUE)), by = quintile]

## Plot (line with points)
p <- ggplot(dom_summ, aes(x = quintile, y = mean_dom, group = 1)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  labs(title = "Mean dominance contribution across quintiles",
       x = "Quintile", y = "Mean dominance") +
  theme_minimal()

print(p)

## (Optional) save
ggsave("/data2/morgante_lab/ukbiobank_projects/prs_effects/results/1000causal_V2/add1000/rep1/pheno/mean_dominance_by_quintile.png",
       p, width = 7, height = 5, dpi = 300)
