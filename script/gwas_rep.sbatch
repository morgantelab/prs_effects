#!/bin/bash
#
#SBATCH --job-name=gwas_rep
#SBATCH --cpus-per-task=20  # Match --threads in PLINK
#SBATCH --partition=fm-bigmem-2
#SBATCH --time=6:00:00
#SBATCH --mem=300G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/gwas_rep_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/gwas_rep_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu



# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=20

# --- Initialize ---

echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="

Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/original_simulation/gwas_rep.R


module unload R/4.2.3

echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="