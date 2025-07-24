#!/bin/bash
#
#SBATCH --job-name=test_prediction_replication
#SBATCH --partition=fm-bigmem-2
#SBATCH --cpus-per-task=40
#SBATCH --mem=800G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/test_prediction_replication_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/test_prediction_replication_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=40

# --- Initialize ---
START_TIME=$(date +%s)
echo "=== Job started at $(date +'%Y-%m-%d %H:%M:%S') ==="

Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/original_simulation/test_prediction.R


module unload R/4.2.3
