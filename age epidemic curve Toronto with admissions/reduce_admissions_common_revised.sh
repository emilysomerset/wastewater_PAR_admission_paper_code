#!/bin/bash
#SBATCH --job-name=reduce_adm
#SBATCH --output=logs/reduce_adm_%j.out
#SBATCH --error=logs/reduce_adm_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load r/4.5.0

cd ~/wastewater_admissions_common_prior_sensitivity

Rscript reduce_admissions_common_revised.R
