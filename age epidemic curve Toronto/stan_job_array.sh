#!/bin/bash
#SBATCH --job-name=stan_wastewater
#SBATCH --account= CENSORED
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=stan_%A_%a.out
#SBATCH --error=stan_%A_%a.err
#SBATCH --array=1-500   # number of IDs in waste_ids.txt

# Load modules
module load r/4.5.0
module load gcc/12.3

# Set CmdStan path for R
export CMDSTAN=/home/somerse5/.cmdstan/cmdstan-2.37.0

# Directory to save results
RESULT_DIR=~/def-justins/wastewater_results
mkdir -p $RESULT_DIR

# Read IDs from file into a Bash array
WASTE_IDS=($(<waste_ids.txt))

# Pick the ID for this array task
ID=${WASTE_IDS[$SLURM_ARRAY_TASK_ID-1]}

# Run the R script
Rscript run_stan_job.R $ID $RESULT_DIR
