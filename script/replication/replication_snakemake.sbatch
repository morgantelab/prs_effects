#!/bin/bash
#
#SBATCH --job-name=prs_effects
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_p_10E-3%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_p_10E-3_%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=2G            
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --partition=fm-bigmem-3
#SBATCH --mail-type=FAIL

# Debugging information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# --- Load required modules ---

echo "Loading modules..."

source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh


# --- Activate Snakemake environment ---
echo "Activating conda environment..."
mamba activate snakemake-7.32

snakemake -p -n -s PRS_scan_Snakefile \
          --configfile PRS_scan_Snakefile \
          --profile slurm \
          --latency-wait 120

# --- Deactivate environment ---    
mamba deactivate


echo "Job finished at: $(date)"