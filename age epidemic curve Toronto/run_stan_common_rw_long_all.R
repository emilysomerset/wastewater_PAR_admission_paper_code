#!/usr/bin/env Rscript
# Longer rerun of the original common-RW case-only model, all 500 wastewater IDs.
suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(magrittr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Expected one SLURM array index")
array_index <- as.integer(args[1])
if (is.na(array_index) || array_index < 1L || array_index > 500L)
  stop("Invalid array index")

project_dir <- "/home/somerse5/wastewater_model"
result_dir <- "/scratch/somerse5/wastewater_model/results_noadm_common_rw_long"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set_cmdstan_path(Sys.getenv("CMDSTAN"))

# Same wastewater draws, order, and seeds as original common-RW run.
set.seed(2917)
ids <- sample(1:3000, 500)
waste_id <- ids[array_index]
seed_i <- 12345L + array_index

# Do not overwrite completed longer fits on an accidental resubmission.
file_name <- file.path(result_dir,
  sprintf("stan_samples_common_rw_long_wastewatersample_%d.RData", waste_id))
if (file.exists(file_name)) {
  message("Already exists; skipping: ", file_name)
  quit(save = "no", status = 0)
}

e <- new.env(parent = globalenv())
load(file.path(project_dir, "epidemic_model_data_noadm.RData"), envir = e)
init <- with(e, case_age_merged$case[
  case_age_merged$earliest_week_end_date == (date_range[1] - 7)])
ratio <- with(e, data_foranalysis_full$analysis_d[[waste_id]] %>%
  arrange(-desc(earliest_week_end_date)) %>%
  filter(earliest_week_end_date >= date_range[1] &
         earliest_week_end_date <= date_range[2]) %>%
  filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed)
stan_data <- with(e, list(
  J = nrow(case_age_matrix),
  I = ncol(case_age_matrix) - 1,
  Y = t(as.matrix(case_age_matrix[, -1])),
  ratio = ratio[seq_len(nrow(case_age_matrix))],
  mean_X0 = init * 2,
  lower_X0 = init,
  sd_X0 = c(1537, 734, 311, 117),
  phi_p = -log(0.5) / 0.4
))
if (length(ratio) < stan_data$J) stop("Not enough wastewater ratio values")
mod <- cmdstan_model(
  stan_file = file.path(project_dir, "model_RWpi_common_noadm.stan"),
  exe_file = file.path(project_dir, "model_RWpi_common_noadm"),
  compile = FALSE
)
fit <- mod$sample(
  data = stan_data,
  seed = seed_i,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  output_dir = result_dir
)
draws <- fit$draws()
wb_rhat <- fit$summary()
sampler_diagnostics <- fit$sampler_diagnostics()
tmp_name <- tempfile(pattern = ".saving_", tmpdir = result_dir, fileext = ".RData")
save(draws, wb_rhat, sampler_diagnostics, file = tmp_name)
if (!file.rename(tmp_name, file_name)) stop("Could not save ", file_name)
unlink(fit$output_files())
message("Saved ", file_name,
        "\nMax Rhat: ", max(wb_rhat$rhat, na.rm = TRUE),
        "\nMinimum bulk ESS: ", min(wb_rhat$ess_bulk, na.rm = TRUE),
        "\nMinimum tail ESS: ", min(wb_rhat$ess_tail, na.rm = TRUE),
        "\nDivergences: ", sum(sampler_diagnostics[, , "divergent__"]))
