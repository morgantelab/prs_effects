#!/bin/bash
#
#SBATCH --job-name=effect_size_estimation_replication
#SBATCH --partition=fm-bigmem-3
#SBATCH --cpus-per-task=54
#SBATCH --mem=800G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/effect_size_estimation_replication_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/effect_size_estimation_replication_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=54

# --- Initialize ---
START_TIME=$(date +%s)
echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="

Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/original_simulation/estimate_effects_allchromosome_combined.R


module unload R/4.2.3