#!/bin/bash
#SBATCH --job-name=parallel_merge       # Job name
#SBATCH --output=parallel_merge_%j.log  # Standard output and error log
#SBATCH --error=parallel_merge_%j.err
#SBATCH --ntasks=1                      # Number of tasks (1 for R script)
#SBATCH --cpus-per-task=1               # Number of cores for parallel processing
#SBATCH --mem=192G                       # Memory per node
#SBATCH --time=24:00:00                 # Max runtime hh:mm:ss

# Load R module
module load r/4.5.0   # adjust based on your cluster R version

# Run R script
Rscript process_within.R

