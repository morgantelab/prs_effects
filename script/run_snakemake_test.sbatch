#!/bin/bash
#
#SBATCH --job-name=prs_effects
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/snakemake/snakemake_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/snakemake/snakemake_%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=10G            
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=your_account
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --partition=fm-bigmem-2
#SBATCH --mail-type=FAIL

# Debugging information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# --- Load required modules ---

echo "Loading modules..."
source /opt/intel/oneapi/mkl/2023.0.0/env/vars.sh intel64
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

module load R/4.2.3
# --- Set project directories ---
PROJECT_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects"
SNAKEMAKE_DIR="$PROJECT_DIR/code/Snakemake_test"
export R_LIB1="/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.2"
export R_LIB2="/opt/ohpc/pub/Software/R-4.2.3/lib64/R/library"

# --- Set thread limits ---
export MKL_NUM_THREADS=4

# Create log directories
mkdir -p /data2/morgante_lab/ukbiobank_projects/prs_effects/results/chr_effects_test

# --- Activate Snakemake environment ---
echo "Activating conda environment..."
mamba activate snakemake-7.32
# Run Snakemake with cluster configuration
PROJECT_DIR="/data2/morgante_lab/ukbiobank_projects/prs_effects"
SNAKEMAKE_DIR="$PROJECT_DIR/code/Snakemake_test"
SNAKEFILE="$SNAKEMAKE_DIR/Snakefile"
cd "$SNAKEMAKE_DIR"
snakemake -p \
    --snakefile "$SNAKEFILE" \
    --cluster-config "config.yaml" \
    --jobs 22 \
    --cluster "sbatch --export=R_LIB1,R_LIB2,ALL"

# --- Deactivate environment ---    
mamba deactivate
module unload R/4.2.3 

echo "Job finished at: $(date)"