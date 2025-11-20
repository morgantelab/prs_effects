#!/bin/bash
#
#SBATCH --job-name=BCT
#SBATCH --output=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/Blood_cell_traits/Blood_cell_traitsp_%j.out
#SBATCH --error=/data2/morgante_lab/ukbiobank_projects/prs_effects/logs/Blood_cell_traits/Blood_cell_traits_%j.err
#SBATCH --time=7-00:00:00       
#SBATCH --mem=200G 
#SBATCH --partition=fm-bigmem-4            
#SBATCH --mail-user=yuxuans@clemson.edu
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=24


/data2/morgante_lab/yuxuans/software/plink2 \
  --bfile /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/filtered/ukb_chr7_filtered \
  --pheno /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_clean.txt \
  --pheno-name WBC_count \
  --covar /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates_plink.txt \
  --covar-col-nums 3-140 \
  --no-input-missing-phenotype \
  --vif 100 \
  --glm hide-covar omit-ref allow-no-covars \
  --threads 24 \
  --out test_WBC

/data2/morgante_lab/yuxuans/software/plink2 \
  --bfile /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/filtered/ukb_chr7_filtered \
  --pheno /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_16traits_clean.txt \
  --pheno-name Monocyte_perc \
  --covar /data2/morgante_lab/ukbiobank_projects/prs_effects/results/Blood_cell_traits/ukb_bloodcell_covariates_plink.txt \
  --covar-col-nums 3-140 \
  --no-input-missing-phenotype \
  --vif 100 \
  --glm hide-covar omit-ref allow-no-covars \
  --threads 24 \
  --out test_MONO

