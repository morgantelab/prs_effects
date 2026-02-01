#!/bin/bash
#
#SBATCH --job-name=BCT_reference
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/BCT_h2/heritibility_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/BCT_h2/heritibility_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G
#SBATCH --partition=fm-bigmem-4
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --mail-type=FAIL

cd /data2/morgante_lab/yuxuans/software/d-ldsc

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh

python get_ldscores.py \
  --bfile /data2/morgante_lab/yuxuans/reference/1kg/all_hg38_hm3_plink1 \
  --out   /data2/morgante_lab/yuxuans/reference/1kg/ldscores_hm3 \
  --ld-wind-cm 1 \
  --yes-really

conda deactivate