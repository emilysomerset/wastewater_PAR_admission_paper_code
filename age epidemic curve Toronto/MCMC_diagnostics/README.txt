Case-only epidemic model, common reporting random walk: MCMC diagnostics

The two CSV files are the finalized parameter-level summaries across 500
conditional fits, based on 4 chains x 2,000 retained draws per fit.

R-hat: median, 75th percentile, maximum, proportions >1.01 and >1.05.
ESS: minimum, 25th percentile, median, separately for bulk and tail ESS.
N is the number of parameter-fit combinations, not the number of fits.

To regenerate on Fir:
  module load r/4.5.0
  Rscript mcmc_diagnostics_case_only.R \
    /scratch/somerse5/wastewater_model/results_noadm_common_rw_long \
    ./MCMC_diagnostics

This requires the original fit RData files (which are not included in GitHub).
The script reads one fit at a time, but the summary table accumulates diagnostics.
The included CSV files are derived from the numerical results supplied in the
research conversation, rounded to the precision shown there.

The case-only divergence count was not supplied for the longer 500-fit run;
no claim about divergences is made in these files.
