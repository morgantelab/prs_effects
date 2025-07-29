#!/bin/bash
#
#SBATCH --job-name=replication_p_10E-3
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_p/replication_p_10E-3%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/replication_p/replication_p_10E-3%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=2G            
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --partition=fm-bigmem-3
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
cd /data2/morgante_lab/ukbiobank_projects/prs_effects/script/replication

snakemake -p \
  -s PRS_scan_Snakefile \
  --configfile PRS_scan_config.yaml \
  --profile /data2/morgante_lab/ukbiobank_projects/prs_effects/script/replication/profiles/slurm \
  --latency-wait 120 \
  --rerun-incomplete 
# --unlock

# --- Deactivate environment ---    
mamba deactivate


echo "Job finished at: $(date)"