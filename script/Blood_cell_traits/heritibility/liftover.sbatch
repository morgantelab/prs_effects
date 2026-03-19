#!/bin/bash
#
#SBATCH --job-name=UKB_liftover
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/Blood_cell_traits/liftover/liftover_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/Blood_cell_traits/liftover/liftover_%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=50G 
#SBATCH --partition=fm-bigmem-4            
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=24


set -euo pipefail

echo "Job started at $(date)"
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
conda activate liftover
echo "Activate conda environment"
LIFTOVER=liftOver

PLINK2=/data2/morgante_lab/yuxuans/software/plink2
CHAIN=/data2/morgante_lab/yuxuans/reference/liftover/hg19ToHg38.over.chain

BED_BIM_DIR=/data2/morgante_lab/data/ukbiobank/genotypes/array
FAM_DIR=/data2/morgante_lab/data/ukbiobank/genotypes/array/fam_04042024
OUT_DIR=/data2/morgante_lab/ukbiobank_projects/prs_effects/results/hg38
mkdir -p ${OUT_DIR}

for CHR in $(seq 1 22); do

  echo "=============================="
  echo "Processing chr${CHR}"
  echo "=============================="

  IN_PREFIX=${BED_BIM_DIR}/ukb22418_c${CHR}_b0_v2_s488243
  FAM=${FAM_DIR}/ukb22418_c${CHR}_b0_v2_s487945.fam
  OUT_PREFIX=${OUT_DIR}/ukb22418_c${CHR}_b0_v2_s488243_hg38

  TMPDIR=liftover_tmp_chr${CHR}
  mkdir -p ${TMPDIR}

  # -------------------------
  # 1) BIM → BED
  # -------------------------
  awk 'BEGIN{OFS="\t"} {print "chr"$1, $4-1, $4, $2}' \
    ${IN_PREFIX}.bim > ${TMPDIR}/chr${CHR}.hg19.bed

  # -------------------------
  # 2) liftOver
  # -------------------------
  ${LIFTOVER} \
    ${TMPDIR}/chr${CHR}.hg19.bed \
    ${CHAIN} \
    ${TMPDIR}/chr${CHR}.hg38.bed \
    ${TMPDIR}/chr${CHR}.unmapped.bed

  # -------------------------
  # 3) update-map
  # -------------------------
  awk 'BEGIN{OFS="\t"} {print $4, $3}' \
    ${TMPDIR}/chr${CHR}.hg38.bed > ${TMPDIR}/update_map_chr${CHR}.txt

    awk '{print $4}' ${TMPDIR}/chr${CHR}.unmapped.bed > ${TMPDIR}/chr${CHR}.unmapped.snplist


  # -------------------------
  # 4) PLINK2 update
  # -------------------------
${PLINK2} \
  --bed ${IN_PREFIX}.bed \
  --bim ${IN_PREFIX}.bim \
  --fam ${FAM} \
  --update-map ${TMPDIR}/update_map_chr${CHR}.txt 1 2 \
  --exclude ${TMPDIR}/chr${CHR}.unmapped.snplist \
  --make-bed \
  --out ${OUT_PREFIX}

  echo "Finished chr${CHR} at $(date)"
  
  rm -rf ${TMPDIR}
done

condat deactivate

echo "All chromosomes finished at $(date)"