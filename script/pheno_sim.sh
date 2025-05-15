#!/bin/bash
#
#SBATCH --job-name=pheno_sim
#SBATCH --partition=fm-bigmem-2
#SBATCH --cpus-per-task=3
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/pheno_sim_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/pheno_sim_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

# --- Load required modules ---
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
module load R/4.2.3

# --- Set MKL threads ---
export MKL_NUM_THREADS=3

# --- Run R script ---u
Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/pheno_sim.R
#Rscript /data2/morgante_lab/ukbiobank_projects/prs_effects/code/pheno_sim_bigsnpr.R
