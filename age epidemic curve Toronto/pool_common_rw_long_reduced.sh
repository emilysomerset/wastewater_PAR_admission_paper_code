#!/bin/bash
#SBATCH --job-name=pool_common_reduced
#SBATCH --output=logs/pool_common_reduced_%j.out
#SBATCH --error=logs/pool_common_reduced_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
set -e
module load r/4.5.0
cd /home/somerse5/wastewater_model
WEEKS_PER_BATCH=5 MAX_FILES=500 Rscript process_noadm_common_rw_long_reduced.R
