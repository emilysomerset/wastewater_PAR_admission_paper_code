#!/bin/bash
#SBATCH --job-name=compile_adm_prior
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=logs/compile_adm_prior_%j.out
#SBATCH --error=logs/compile_adm_prior_%j.err
set -euo pipefail
module load r/4.5.0
module load gcc/12.3
export CMDSTAN=/home/somerse5/.cmdstan/cmdstan-2.37.0
cd /home/somerse5/wastewater_admissions_common_prior_sensitivity
Rscript -e 'library(cmdstanr); set_cmdstan_path(Sys.getenv("CMDSTAN")); cmdstan_model("model_admissions_common_revised_priors.stan")'
