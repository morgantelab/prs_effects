#!/bin/bash
#
#SBATCH --job-name=replication_p
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_centered/replication_p_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_centered/replication_p_%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=2G 
#SBATCH --partition=fm-bigmem-4            
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --mail-type=FAIL

# Debugging information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"

# --- Load required modules ---
echo "Loading modules..."
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/conda.sh
source /opt/ohpc/pub/Software/mamba-rocky/etc/profile.d/mamba.sh

export R_LIB1="/data2/morgante_lab/yuxuans/software/R/x86_64-pc-linux-gnu-library/4.2"
export R_LIB2="/opt/ohpc/pub/Software/R-4.2.3/lib64/R/library"

# --- Activate Snakemake environment ---
echo "Activating conda environment..."
mamba activate snakemake-7.32
cd /data2/morgante_lab/ukbiobank_projects/prs_effects/script/replication_centered

# Unlock first (safe to run even if not locked)
snakemake \
  -s PRS_snakefile \
  --configfile PRS_scan_config.yaml \
  --profile profiles/slurm \
  --unlock

snakemake \
  -p \
  -s PRS_snakefile \
  --configfile PRS_scan_config.yaml \
  --profile profiles/slurm \
  --latency-wait 120


# --- Deactivate environment ---    
mamba deactivate

echo "Job finished at: $(date)"