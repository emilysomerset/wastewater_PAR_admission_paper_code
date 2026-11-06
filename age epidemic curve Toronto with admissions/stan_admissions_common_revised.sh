#!/bin/bash
#SBATCH --job-name=stan_adm_prior
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-500
#SBATCH --output=logs/stan_adm_prior_%A_%a.out
#SBATCH --error=logs/stan_adm_prior_%A_%a.err
set -euo pipefail
module load r/4.5.0
module load gcc/12.3
export CMDSTAN=/home/somerse5/.cmdstan/cmdstan-2.37.0
export OMP_NUM_THREADS=1
cd /home/somerse5/wastewater_admissions_common_prior_sensitivity
Rscript run_admissions_common_revised.R "$SLURM_ARRAY_TASK_ID"
