# Age-Structured SARS-CoV-2 Epidemic Reconstruction in Toronto

This repository contains the R and Stan code used to reconstruct age-specific SARS-CoV-2 infections in Toronto using wastewater surveillance data, reported COVID-19 cases, and hospital admissions.

Two Bayesian age-structured epidemic models are considered:

1. **Case-only model:** Reconstructs latent infections using reported COVID-19 cases and wastewater-derived reproduction rates.
2. **Cases and hospital admissions model:** Extends the case-only model by incorporating age-specific hospital admissions to provide additional information about latent infections.

Both models use a modular Bayesian inference approach in which wastewater-derived reproduction rates are treated as fixed within each conditional epidemic model fit. Uncertainty from the wastewater model is propagated by fitting the epidemic models separately to 500 posterior samples of the wastewater-derived reproduction rates.

## Model descriptions

### Case-only model

The case-only model reconstructs weekly latent infections for four age groups: 0–39, 40–59, 60–79, and 80 years and older.

The model incorporates age mixing and a common time-varying reporting probability, modelled as a Gaussian random walk on the logit scale.

The corresponding code and results are contained in `age epidemic curve Toronto/`.

### Cases and hospital admissions model

The second model extends the case-only framework by incorporating weekly age-specific hospital admissions.

Hospitalization probabilities vary over time and between age groups. Their temporal dependence is modelled using age-specific AR(1) processes on the logit scale, with a shared autoregressive correlation and innovation standard deviation.

The model also includes a common time-varying reporting probability for non-hospitalized infections.

The observation model distinguishes three mutually exclusive outcomes: hospitalized infections, reported non-hospitalized infections, and unreported non-hospitalized infections.

The corresponding code and results are contained in `age epidemic curve Toronto with admissions/`.

## Repository organization

### `age epidemic curve Toronto/`

Contains the case-only model, including:

- `1_make_epidemic_data.R` — Prepares epidemic model inputs.
- `epidemic_model_data_noadm.RData` — Prepared epidemic model data.
- `model_RWpi_common_noadm.stan` — Stan model specification.
- `run_stan_common_rw_long_all.R` — Fits the conditional epidemic models.
- `stan_common_rw_long_all.sh` — SLURM array submission script.
- `reduce_common_rw_long.R` and `reduce_common_rw_long.sh` — Reduce the number of posterior draws used for pooling.
- `process_noadm_common_rw_long_reduced.R` — Calculates pooled posterior summaries.
- `pool_common_rw_long_reduced.sh` — Submits the pooling job.
- `results_nolag_toronto_pooled_common_rw.RData` — Final pooled posterior summaries.
- `MCMC_diagnostics/` — Convergence diagnostic code and results.
- `prior_posterior/` — Posterior draws and code for prior–posterior comparisons.
- `plots_for_manuscript.R` — Generates figures used in the manuscript.
- `plots/` — Generated figures.

### `age epidemic curve Toronto with admissions/`

Contains the cases and hospital admissions model, including:

- `1_make_epidemic_data.R` — Prepares epidemic model inputs.
- `epidemic_model_nonlag_data.RData` — Prepared epidemic model data.
- `model_admissions_common_revised_priors.stan` — Stan model specification.
- `run_admissions_common_revised.R` — Fits the conditional epidemic models.
- `compile_admissions_common_revised.sh` — Compiles the Stan model.
- `stan_admissions_common_revised.sh` — SLURM array submission script.
- `reduce_admissions_common_revised.R` and `reduce_admissions_common_revised.sh` — Reduce posterior draws for pooling.
- `process_admissions_common_revised_reduced.R` — Calculates pooled posterior summaries.
- `pool_admissions_revised_reduced.sh` — Submits the pooling job.
- `admissions_common_revised_pooled_reduced.RData` — Final pooled posterior summaries.
- `MCMC_diagnostics/` — Convergence diagnostic code and results.
- `prior_posterior/` — Posterior draws and code for prior–posterior comparisons.
- `plots_for_manuscript.R` — Generates figures used in the manuscript.
- `plots/` — Generated figures.

### `wastewater_model`

Contains the code to fit the wastepaper model including:
- `1_wastewatermodel_CAN.R` - Fits the wastewater model to all of Canada (we filter for the results from Toronto for paper)
- `wastewater_model_CAN_phachist.RData` - Fitted wastewater model
- `2_plots_forpaper.R` - Generates figures used in manuscript.

Will call on code from
- `cpp/` - C++ template for wastewater model
- `data/` - Data downloaded from various sources
- `functions_general/` - General helper functions

## Model fitting and posterior inference

Both epidemic models were fitted using Hamiltonian Monte Carlo (HMC) through Stan.

For each model, 500 conditional fits were performed, corresponding to 500 posterior samples from the wastewater model.

Each conditional fit used:

| Setting | Value |
|---|---|
| MCMC chains | 4 |
| Warmup iterations per chain | 1,000 |
| Sampling iterations per chain | 2,000 |
| Retained posterior draws per fit | 8,000 |
| Wastewater posterior samples | 500 |

To reduce computational requirements during posterior processing, 250 draws were selected from each chain, yielding 1,000 draws per conditional fit.

The resulting pooled posterior summaries were calculated from 500,000 draws per age-week.

Posterior summaries include medians and credible intervals for latent infections, cumulative infections, reporting probabilities, and, for the admissions model, hospitalization probabilities.

Cumulative infections were calculated within each posterior draw before computing posterior quantiles.

## Prior–posterior comparisons

The `prior_posterior/` directories contain the data and R scripts used to compare prior and posterior distributions of selected parameters.

For the case-only model, the comparison focuses on the reporting random-walk innovation standard deviation.

For the admissions model, comparisons include the reporting and hospitalization innovation standard deviations, hospitalization AR(1) correlation, and age-specific baseline hospitalization probabilities.

The figures use analytical prior densities and kernel density estimates of posterior samples. Posterior draws are combined across the 500 conditional fits.

## Reproducing the analysis

The general workflow for each model is:

1. Prepare the epidemic model data and obtain the wastewater posterior samples.
2. Compile the Stan model.
3. Fit the 500 conditional epidemic models using the SLURM array scripts.
4. Calculate MCMC convergence diagnostics.
5. Reduce the posterior draws to 1,000 per conditional fit.
6. Pool the posterior draws and calculate posterior summaries.
7. Generate prior–posterior comparisons and manuscript figures.

The fitting and posterior processing were performed on the Digital Research Alliance of Canada's Fir computing cluster.

## Software

The final epidemic models were fitted using R 4.5.0 and CmdStan 2.37.0 through the `cmdstanr` interface.

R packages used throughout the analysis include `dplyr`, `ggplot2`, and `tidyr`.

## Reproducibility and data availability

The repository contains the Stan model specifications, fitting and pooling scripts, convergence diagnostic code, and final pooled posterior summaries.

The complete posterior draws from all 500 conditional fits are not included because of their size.

Some scripts contain paths specific to the Fir computing environment and must be adapted when running the analysis on another system.