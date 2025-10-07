#!/usr/bin/env Rscript
# PRS parity check: manual g*beta (NA->0) vs PLINK2 --score (cols=+scoresums, no-mean-imputation)
# - Uses clumped index SNPs + GWAS ADD betas
# - Extracts TEST samples only
# - Compares per-chromosome and genome-wide parity

suppressPackageStartupMessages({library(data.table); library(genio)})

## =========================
## Config (edit paths if needed)
## =========================
rep <- 2
root <- "/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication_scenario4"
pheno_dir <- file.path(root, "pheno")
gwas_dir  <- file.path(root, "gwas")
out_dir   <- file.path(root, "test")
geno_pref <- file.path(pheno_dir, "ukb_chr")   # => ukb_chr{CHR}_filtered
plink2    <- "/data2/morgante_lab/yuxuans/software/plink2"

test_ids_file  <- file.path(pheno_dir, sprintf("test_ID_rep%d.txt", rep))
combined_clumps<- file.path(out_dir, "clumping/0.01",
                            sprintf("gwas_clump_rep%d_r2_0.5.allchr.clumps", rep))
gwas_all       <- file.path(gwas_dir, sprintf("gwas_rep%d.allchr.glm.linear", rep))

tmp <- file.path(out_dir, sprintf("prs_parity_rep%d_tmp", rep)); dir.create(tmp, TRUE)
parity_out <- file.path(out_dir, sprintf("prs_parity_rep%d.tsv", rep))

## =========================
## Helpers
## =========================
norm_id <- function(x) as.character(trimws(x))
get_bim_alleles <- function(bim){
  nm <- names(bim)
  if (all(c("alt","ref") %in% nm)) return(list(a1=toupper(bim$alt), a2=toupper(bim$ref)))
  if (all(c("a1","a2") %in% nm))   return(list(a1=toupper(bim$a1),  a2=toupper(bim$a2)))
  if (all(c("A1","A2") %in% nm))   return(list(a1=toupper(bim$A1),  a2=toupper(bim$A2)))
  stop("BIM missing allele columns (need alt/ref or a1/a2 or A1/A2)")
}
is_comp <- function(a,b){ comp<-c(A="T",T="A",C="G",G="C"); !is.na(a)&&!is.na(b)&&!is.na(comp[a])&&comp[a]==b }

## =========================
## Load inputs
## =========================
# TEST IIDs (2-col keep => IID=col2; 1-col => col1)
tid <- fread(test_ids_file, header=FALSE)
test_iids <- norm_id(if (ncol(tid)>=2) tid[[2]] else tid[[1]])
cat("TEST IIDs:", length(test_iids), "\n")

# clumped index SNPs (CHR, SNP)
cl <- fread(combined_clumps)
if ("#CHROM" %in% names(cl)) setnames(cl, "#CHROM", "CHR")
if ("ID" %in% names(cl) & !"SNP" %in% names(cl)) setnames(cl, "ID", "SNP")
stopifnot(all(c("CHR","SNP") %in% names(cl)))
cl[, CHR:=as.integer(CHR)][, SNP:=as.character(SNP)]
cl <- unique(cl[,.(CHR,SNP)])
cat("Index SNPs in .clumps:", nrow(cl), "\n")

# GWAS ADD betas (CHR, SNP, A1, BETA)
gw <- if (file.exists(gwas_all)) fread(gwas_all) else {
  rbindlist(lapply(sprintf("%s/gwas_rep%d_chr%d.pheno.glm.linear", gwas_dir, rep, 1:22),
                   function(p) if (file.exists(p)) fread(p) else NULL),
            fill=TRUE)
}
if ("#CHROM" %in% names(gw)) setnames(gw,"#CHROM","CHR")
if ("ID" %in% names(gw) & !"SNP" %in% names(gw)) setnames(gw,"ID","SNP")
stopifnot(all(c("CHR","SNP","A1","TEST","BETA") %in% names(gw)))
gw[, CHR:=as.integer(CHR)][, SNP:=as.character(SNP)]
gw <- gw[TEST=="ADD", .(CHR,SNP,A1,BETA)]
gw <- unique(gw, by=c("CHR","SNP"))
cat("GWAS ADD rows:", nrow(gw), "\n")

# Intersection: clumped index SNPs with ADD betas
idx_beta <- merge(cl, gw, by=c("CHR","SNP"))
cat("Index SNPs with ADD betas:", nrow(idx_beta), "\n")

## =========================
## Accumulators
## =========================
prs_hand  <- numeric(length(test_iids)); names(prs_hand)  <- test_iids
prs_plink <- numeric(length(test_iids)); names(prs_plink) <- test_iids

## =========================
## Main loop per chromosome
## =========================
for (chr in 1:22) {
  sb <- idx_beta[CHR==chr]
  if (!nrow(sb)) { cat(sprintf("chr%02d: no SNPs\n", chr)); next }
  
  # Write score file (SNP A1 BETA)
  scoref  <- file.path(tmp, sprintf("chr%d.score", chr))
  fwrite(sb[,.(SNP,A1,BETA)], scoref, sep="\t")
  
  # Extract samples + SNPs
  bfile   <- sprintf("%s%d_filtered", geno_pref, chr)
  outp    <- file.path(tmp, sprintf("chr%d_extract", chr))
  snplist <- file.path(tmp, sprintf("chr%d.snplist", chr)); fwrite(sb[,.(SNP)], snplist, col.names=FALSE)
  
  system2(plink2, c("--bfile", bfile, "--keep", test_ids_file, "--extract", snplist,
                    "--make-bed", "--out", outp, "--threads","8"))
  
  if (!file.exists(paste0(outp,".bed"))) { cat(sprintf("chr%02d: extract failed\n", chr)); next }
  
  ## --- PLINK --score (SUM, no-mean-imputation) ---
  scr_out <- file.path(tmp, sprintf("chr%d_score", chr))
  system2(plink2, c("--bfile", outp,
                    "--score", scoref, "1","2","3", "no-mean-imputation", "cols=+scoresums",
                    "--out", scr_out, "--threads","4"))
  
  ssc <- fread(paste0(scr_out, ".sscore"))
  sumcol <- grep("^SCORE\\d*_SUM$|^SCORESUM$|_SUM$", names(ssc), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(sumcol)) stop("No *_SUM column in .sscore; upgrade plink2 or add cols=+scoresums")
  ssc[, IID:=norm_id(IID)]; setkey(ssc, IID)
  prs_plink_chr <- ssc[J(test_iids)][[sumcol]]
  prs_plink[!is.na(prs_plink_chr)] <- prs_plink[!is.na(prs_plink_chr)] + prs_plink_chr[!is.na(prs_plink_chr)]
  
  ## --- Manual g * beta (NA -> 0; align to GWAS A1, flip when A1==BIM A2/comp) ---
  gp  <- genio::read_plink(outp, verbose=FALSE)
  fam <- norm_id(gp$fam$id)
  X   <- t(gp$X)                           # indiv x snp
  ord <- match(test_iids, fam); if (any(is.na(ord))) stop("FAM vs TEST mismatch")
  X <- X[ord, , drop=FALSE]
  
  bim  <- as.data.table(gp$bim)
  alle <- get_bim_alleles(bim)
  
  # map GWAS rows to BIM order
  idx <- match(bim$id, sb$SNP)
  A1  <- toupper(sb$A1[idx])
  B   <- sb$BETA[idx]
  
  # keep only SNPs with valid allele mapping (A1 matches BIM A1/A2 or complements)
  a <- alle$a1; b <- alle$a2
  ok <- !is.na(A1) & !is.na(B) &
    (A1==a | A1==b | mapply(is_comp, A1, a) | mapply(is_comp, A1, b))
  if (!any(ok)) { cat(sprintf("chr%02d: all SNPs dropped by allele check\n", chr)); next }
  
  # flip when GWAS A1 corresponds to BIM A2 (or its complement)
  flip <- (A1==b) | mapply(is_comp, A1, b)
  
  # subset to ok columns
  X   <- X[, ok, drop=FALSE]
  B   <- B[ok]
  flip<- flip[ok]
  
  # NA -> 0 to mirror no-mean-imputation
  X[is.na(X)] <- 0
  if (any(flip)) X[, which(flip)] <- 2 - X[, which(flip), drop=FALSE]
  
  contrib_hand <- as.numeric(X %*% B)
  prs_hand <- prs_hand + contrib_hand
  
  # optional per-chr parity
  ssc_chr <- ssc[J(test_iids)][[sumcol]]
  cat(sprintf("chr%02d: SNPs_kept=%d, flips=%d, corr(hand_chr, plink_chr)=%.6f\n",
              chr, ncol(X), sum(flip),
              suppressWarnings(cor(contrib_hand, ssc_chr, use="complete.obs"))))
}

## =========================
## Compare genome-wide & save
## =========================
okv <- !(is.na(prs_hand) | is.na(prs_plink))
cat("Genome-wide cor(hand, plink) =", suppressWarnings(cor(prs_hand[okv], prs_plink[okv])), "\n")
cat("Max abs diff =", max(abs(prs_hand[okv]-prs_plink[okv])), "\n")

out_tab <- data.table(IID=test_iids, PRS_hand=prs_hand, PRS_plink=prs_plink,
                      diff = prs_hand - prs_plink)
fwrite(out_tab, parity_out, sep="\t")
cat("Wrote:", parity_out, "\n")

# Optional cleanup
unlink(tmp, recursive=TRUE)

cor(contrib_hand, ssc_chr, use="complete.obs")
summary(contrib_hand - ssc_chr)

# 看 SNP 层面对齐情况
cbind(bim$id, bim$alt, bim$ref, A1, flip)[1:20,]
