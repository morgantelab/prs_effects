#!/bin/bash
#
#SBATCH --job-name=prs_snakemake
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-2
#SBATCH --time=12:30:00
#SBATCH --mem=10G
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/snakemake/snakemake_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/errors/snakemake/snakemake_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxuans@clemson.edu

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
SNAKEMAKE_DIR="$PROJECT_DIR/code/Snakemake"
export R_LIB1="/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.2"
export R_LIB2="/opt/ohpc/pub/Software/R-4.2.3/lib64/R/library"

# --- Set thread limits ---
export MKL_NUM_THREADS=4
# --- Create directories ---
mkdir -p "$PROJECT_DIR/temp"
mkdir -p "$PROJECT_DIR/logs/"{process_chromosome,combine_results,generate_plots}

SNAKEFILE="$SNAKEMAKE_DIR/Snakefile"
DIRS_CONFIG="$SNAKEMAKE_DIR/dirs.yaml"
CLUSTER_CONFIG="$SNAKEMAKE_DIR/cluster_config.yaml"
# --- Activate Snakemake environment ---
echo "Activating conda environment..."
mamba activate snakemake-7


# Run Snakemake
cd "$SNAKEMAKE_DIR"
snakemake -p \
    --snakefile "$SNAKEFILE" \
    --cluster-config "config.yaml" \
    --jobs 1 \
    --cluster "sbatch --export=R_LIB1,R_LIB2,ALL"

# --- Deactivate environment ---
mamba deactivate
module unload R/4.2.3 

echo "Job finished at: $(date)"