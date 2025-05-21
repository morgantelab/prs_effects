#!/bin/bash
#
#SBATCH --job-name=plink2_gwas
#SBATCH --cpus-per-task=4  # Match --threads in PLINK
#SBATCH --partition=fm-bigmem-2
#SBATCH --time=6:00:00
#SBATCH --mem=30G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/plink2_GWAS_lare_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/plink2_GWAS_large_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Configuration ---
BED_BIM_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects/filter"
OUTPUT_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects/results"
PLINK2="/data2/morgante_lab/yuxuans/software/plink2"
PHENO_FILE="${OUTPUT_DIR}/simulated_pheno_large.txt"

# --- Initialize timing ---
START_TIME=$(date +%s)
echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="

# --- Create output directory ---
mkdir -p "${OUTPUT_DIR}/gwas_results"

# --- Process chromosomes ---
for CHR in {1..22}; do
    CHR_START=$(date +%s)
    echo "[$(date +'%H:%M:%S')] Processing chromosome ${CHR}"
    
    ${PLINK2} \
        --bfile "${BED_BIM_DIR}/ukb_chr${CHR}_filtered" \
        --pheno "${PHENO_FILE}" \
        --glm 'no-x-sex' 'allow-no-covars'  \
        --out "${OUTPUT_DIR}/gwas_results/chr${CHR}_gwas_large" \
        --threads 4
    
    CHR_DURATION=$(( $(date +%s) - CHR_START ))
    echo "[$(date +'%H:%M:%S')] Completed chromosome ${CHR} in $(printf '%02d:%02d:%02d\n' $((CHR_DURATION/3600)) $((CHR_DURATION%3600/60)) $((CHR_DURATION%60)))"
done


# --- Calculate total runtime ---
END_TIME=$(date +%s)
TOTAL_DURATION=$(( END_TIME - START_TIME ))
echo "=== Job completed at $(date +'%Y-%m-%d %H:%M:%S') ==="
echo "Total runtime: $(printf '%02d:%02d:%02d\n' $((TOTAL_DURATION/3600)) $((TOTAL_DURATION%3600/60)) $((TOTAL_DURATION%60)))"

