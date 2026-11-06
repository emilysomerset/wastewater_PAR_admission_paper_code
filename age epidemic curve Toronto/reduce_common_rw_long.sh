#!/bin/bash
#SBATCH --job-name=reduce_common
#SBATCH --output=logs/reduce_common_%j.out
#SBATCH --error=logs/reduce_common_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
set -e
module load r/4.5.0
cd /home/somerse5/wastewater_model
Rscript reduce_common_rw_long.R
