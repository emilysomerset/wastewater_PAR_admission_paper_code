MCMC diagnostics: revised admissions + cases model, common reporting RW

Run on Fir:
  module load r/4.5.0
  Rscript mcmc_diagnostics_admissions.R

By default the script reads the 500 original fit files from:
  /scratch/somerse5/wastewater_admissions_common_prior_sensitivity/fits
and writes the three CSVs in the parent directory.

To change paths, set ADMISSIONS_FIT_DIR and DIAGNOSTICS_OUTPUT_DIR.

The committed summary CSVs are based on the 500 longer conditional fits
(4 chains x 2,000 retained iterations). N counts parameter-fit combinations.
The proportion columns are fractions, not percentages.
