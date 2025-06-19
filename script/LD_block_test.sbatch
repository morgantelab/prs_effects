#!/bin/bash
#
#SBATCH --job-name=LD_blcok_test
#SBATCH --partition=fm-bigmem-2
#SBATCH --cpus-per-task=20
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/ld_block_test_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/ld_block_test_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=20

# --- Run R script ---
Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/LD_block_test.R
