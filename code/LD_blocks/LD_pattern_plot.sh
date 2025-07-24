#!/bin/bash
#
#SBATCH --job-name=chr1_LD_pattern
#SBATCH --cpus-per-task=24
#SBATCH --partition=fm-bigmem-4
#SBATCH --time=6:00:00
#SBATCH --mem=200G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/chr1_LD_pattern_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/chr1_LD_pattern_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=24

# --- Run R script ---
Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/LD_blocks/LD_patterns.R