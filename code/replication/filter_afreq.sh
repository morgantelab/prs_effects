#!/bin/bash
#
#SBATCH --job-name=plink2_filter_all_chr
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-2
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/plink2_all_chr_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/plink2_all_chr_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# ==============================================
# PLINK 2.0 Multi-Chromosome Filtering Script
# ==============================================

# --- Configuration ---
BED_BIM_DIR="/data2/morgante_lab/data/ukbiobank/genotypes/array"
FAM_DIR="/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024"
OUTPUT_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/filter"
ID_FILE="/data2/morgante_lab/ukbiobank_projects/prs_effects/results/replication/ID_filtered.txt"
PLINK2="/data2/morgante_lab/yuxuans/software/plink2"

# --- Initialize ---
START_TIME=$(date +%s)
echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="
echo "Output will be written to: ${OUTPUT_DIR}"

source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/sample_ini_selection.R


# --- Chromosome Loop ---
for CHR in {1..22}; do
    CHR_START=$(date +%s)
    echo "[$(date +'%H:%M:%S')] Processing chromosome ${CHR}..."
    
    # Input files
    BED="${BED_BIM_DIR}/ukb22418_c${CHR}_b0_v2_s488243.bed"
    BIM="${BED_BIM_DIR}/ukb22418_c${CHR}_b0_v2_s488243.bim"
    FAM="${FAM_DIR}/ukb22418_c${CHR}_b0_v2_s487945.fam"

    # Run PLINK2
    "${PLINK2}" \
      --bed "${BED}" \
      --bim "${BIM}" \
      --fam "${FAM}" \
      --keep "${ID_FILE}" \
      --rm-dup exclude-all \
      --snps-only just-acgt \
      --maf 0.01 \
      --mac 5 \
      --geno 0.1 \
      --hwe 1e-5 \
      --threads 8 \
      --make-bed \
      --out "${OUTPUT_DIR}/ukb_chr${CHR}_filtered"
done

# --- Get Allele Frequency ---
for CHR in {1..22} ; do
    "${PLINK2}" \
      --bfile "${OUTPUT_DIR}/ukb_chr${CHR}_filtered" \
      --freq \
      --out "${OUTPUT_DIR}/ukb_chr${CHR}_filtered"
done


