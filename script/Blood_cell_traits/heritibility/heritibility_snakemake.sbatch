#!/bin/bash
#
#SBATCH --job-name=BCT_h2
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/BCT_h2/heritibility_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/BCT_h2/heritibility_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=8G
#SBATCH --partition=fm-bigmem-4
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --mail-type=FAIL

set -euo pipefail

echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"

# --- Load required conda/mamba init ---
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

# --- Activate Snakemake environment ---
export R_LIB1="/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.2"
export R_LIB2="/opt/ohpc/pub/Software/R-4.2.3/lib64/R/library"

mamba activate snakemake-7.32

###############################################################################
# EDIT HERE: set your workflow directory + file names
###############################################################################
WORKDIR="/data2/morgante_lab/ukbiobank_projects/prs_effects/script/Blood_cell_traits/heritibility"   # <- h2 Snakefile 所在目录
SNAKEFILE="h2_snakefile"                                                               # <- Snakefile 文件名
CONFIG="heritibility.yaml"                                                             # <- config 文件名
PROFILE="profiles/slurm"                                                               # <- slurm profile
###############################################################################

cd "${WORKDIR}"

echo "Running workflow in: $(pwd)"
echo "Snakefile: ${SNAKEFILE}"
echo "Config: ${CONFIG}"

# Unlock (safe even if not locked)
snakemake \
  -s "${SNAKEFILE}" \
  --configfile "${CONFIG}" \
  --profile "${PROFILE}" \
  --unlock

# Run
snakemake \
  -p \
  -s "${SNAKEFILE}" \
  --configfile "${CONFIG}" \
  --profile "${PROFILE}" \
  --latency-wait 120

mamba deactivate
echo "Job finished at: $(date)"
