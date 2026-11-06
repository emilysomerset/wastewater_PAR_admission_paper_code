#!/bin/bash
#SBATCH --job-name=pool_adm_revised
#SBATCH --output=logs/pool_adm_revised_%j.out
#SBATCH --error=logs/pool_adm_revised_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

module load r/4.5.0

cd ~/wastewater_admissions_common_prior_sensitivity

MAX_FILES=500 WEEKS_PER_BATCH=2 \
  Rscript process_admissions_common_revised_reduced.R
