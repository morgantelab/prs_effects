#!/bin/bash
#
#SBATCH --job-name=effect_size_estimation
#SBATCH --partition=fm-bigmem-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=800G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/effect_size_estimation_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/effect_size_estimation_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=4

# --- Run R script ---u
Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/effect_size_estimation.R