#!/bin/bash
#SBATCH --job-name=stan_common_long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/stan_common_long_%A_%a.out
#SBATCH --error=logs/stan_common_long_%A_%a.err
#SBATCH --array=1-500

module load r/4.5.0
module load gcc/12.3
export CMDSTAN=/home/somerse5/.cmdstan/cmdstan-2.37.0
export OMP_NUM_THREADS=1
cd /home/somerse5/wastewater_model || exit 1
Rscript run_stan_common_rw_long_all.R "$SLURM_ARRAY_TASK_ID"
