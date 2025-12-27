#!/usr/bin/env Rscript

# ===========================
# run_stan_job.R
# ===========================
# Runs your wastewater Stan model for a given sample ID
# Arguments:
#   1. waste_id (integer ID)
#   2. RESULT_DIR (where to save results)
# ===========================

# ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Need two arguments: waste_id and RESULT_DIR")
waste_id <- as.integer(args[1])
RESULT_DIR <- args[2]

# ---- Load packages ----
library(cmdstanr)
library(rstan)
library(dplyr)
library(magrittr)
library(reshape2)
library(posterior)
library(lubridate)
library(tidyr)

# ---- Set CmdStan path ----
cmdstan_path <- Sys.getenv("CMDSTAN")
if (cmdstan_path == "") stop("CMDSTAN environment variable is not set!")
set_cmdstan_path(cmdstan_path)
message("Using CmdStan version: ", cmdstan_version())

# ---- Load data ----
load(file.path(RESULT_DIR, "epidemic_model_nonlag_data.RData"))

# ---- Initialize infections ----
init <- case_age_merged$case[case_age_merged$earliest_week_end_date == (date_range[1]-7)]

# ---- Load model ----
mod <- cmdstan_model(file.path(RESULT_DIR, "model_final_switchC_RWpi_ARtau_noforecast_noendemic.stan"))

# ---- Deterministic seed ----
base_seed <- 12345L
seed_i <- base_seed + waste_id

# ---- Prepare stan_data ----
ratio <- data_foranalysis_full$analysis_d[[waste_id]] %>%
  arrange(-desc(earliest_week_end_date)) %>%
  filter(earliest_week_end_date >= date_range[1] & earliest_week_end_date <= date_range[2]) %>%
  filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed

stan_data <- list(
  J = nrow(case_age_matrix),
  I = ncol(case_age_matrix) - 1,
  J_future = length(ratio),
  Y = t(as.matrix(case_age_matrix[,-1])),
  ratio = ratio,
  mean_X0 = init*2,
  lower_X0 = init,
  sd_X0 = c(786,233,124,30),
  phi_p = -log(0.5)/0.4,
  A = t(as.matrix(hosp_age_matrix[,-1])),
  phi_tau = -log(0.5)/0.4,
  tau_mean_prior = c(-4,-2,-2,-2)
)

# ---- Initialize function ----
pi_guess <- 0.5
init_fun <- function() {
  list(
    X_init = pmax(init*2 + rnorm(length(init),0, stan_data$sd_X0), init),
    X = pmax(
      stan_data$Y + rpois(length(stan_data$Y), lambda = stan_data$Y / pi_guess),
      stan_data$Y + 1
    )
  )
}

# ---- Run Stan ----
fit <- mod$sample(
  data = stan_data,
  seed = seed_i,
  chains = 6,
  parallel_chains = 6,
  show_messages = TRUE,
  iter_warmup = 4000,
  iter_sampling = 15000,
  save_warmup = FALSE,
  init = init_fun,
  thin = 30
)

# ---- Extract draws and summary ----
draws <- fit$draws()
wb_rhat <- fit$summary()

# ---- Save results ----
file_name <- file.path(RESULT_DIR, paste0("stan_samples_nolag_wastewatersample_", waste_id, ".rds"))
save(file = file_name, list = c("draws", "wb_rhat"))
message("Saved results to: ", file_name)
