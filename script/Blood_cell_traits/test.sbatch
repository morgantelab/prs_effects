#!/bin/bash
#SBATCH --job-name=kvik_step2_test
#SBATCH --output=kvik_step2_test.out
#SBATCH --error=kvik_step2_test.err
#SBATCH --partition=fm-bigmem-4
#SBATCH --account=yuxuans
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G


echo "==== LDAK KVIK Step1 Test ===="

LDAK="/data2/morgante_lab/yuxuans/software/ldak6.1.linux"
OUT_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits"
MERGED="${OUT_DIR}/merged/ukb_merged"
PHENO="${OUT_DIR}/ukb_bloodcell_16traits_quantnorm.txt"
COVAR="${OUT_DIR}/ukb_bloodcell_covariates.txt"
KEEP="${OUT_DIR}/ukb_bloodcell_final_keep.txt"

mkdir -p ${OUT_DIR}/kvik_test

${LDAK} \
  --kvik-step2 ${OUT_DIR}/kvik_test/kvik_output \
  --bfile ${MERGED} \
  --pheno ${PHENO} \
  --covar ${COVAR} \
  --mpheno ALL \
  --max-threads 32

echo "==== Test Completed ===="
echo "Output directory:"
echo "    ${OUT_DIR}/kvik_test/"
echo "List files:"
ls -lh ${OUT_DIR}/kvik_test/
