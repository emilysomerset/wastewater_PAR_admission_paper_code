# Packages
library(rstan) # rstan_2.32.7 
library(dplyr) # dplyr_1.1.4
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(doParallel) # doParallel_1.0.17
library(cmdstanr) # cmdstanr_0.9.0
library(posterior) #posterior_1.6.0
library(lubridate) # lubridate_1.9.4
library(tidyr)


# Read data
load("epidemic_model_data_noadm.RData")

init = case_age_merged$case[case_age_merged$earliest_week_end_date == (date_range[1]-7)] 


# Read model
mod <- cmdstan_model("model_final_switchC_RWpi_ARtau_noforecast_noendemic_noadm_diffX0prior.stan")

# Make a list of end points for the segments 

n_wastewater <- 500
set.seed(2917)
ss <- sample(1:3000,500)

grid <- expand.grid(
  waste_id = ss,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) 

# Deterministic Stan seeds (one per job)
base_seed <- 12345L
stan_seeds <- base_seed + seq_len(nrow(grid))

# Start cluster

myCluster <- makeCluster(10, # number of cores to use
                         type = "PSOCK")

registerDoParallel(myCluster)
print(getDoParWorkers())

# record start time
start_time <- Sys.time()

res <- foreach(i = seq_len(nrow(grid))) %dopar% {
  # Packages
  library(rstan) # rstan_2.32.7 
  library(dplyr) # dplyr_1.1.4
  library(magrittr) # magrittr_2.0.3
  library(reshape2) # reshape2_1.4.4
  library(doParallel) # doParallel_1.0.17
  library(cmdstanr) # cmdstanr_0.9.0
  library(posterior) #posterior_1.6.0
  library(lubridate) # lubridate_1.9.4
  library(tidyr)

  waste_id <- grid$waste_id[i]
  seed_i <- stan_seeds[i]
  
  ratio = data_foranalysis_full$analysis_d[[waste_id]] %>% 
    arrange(-desc(earliest_week_end_date)) %>% 
    filter(earliest_week_end_date >= date_range[1] & earliest_week_end_date <= date_range[2]) %>% 
    filter(!is.na(ratio_cases)) %$% ratio_v_u_fixed

  
  stan_data <- list(
    J = nrow(case_age_matrix),
    I = ncol(case_age_matrix)-1,
    J_future = length(ratio),
    
    # Infection model  
    Y = as.matrix(case_age_matrix[,-1]) %>% t(),
    ratio = ratio,
    mean_X0 = init*2,
    lower_X0 = init,
    sd_X0 = c(1537, 734, 311, 117),
    phi_p = -log(0.5)/0.4
  ) 
  
  pi_guess <- 0.5
  init_fun <- function() {
    list(
      X = pmax(
        stan_data$Y + rpois(length(stan_data$Y), lambda = stan_data$Y /pi_guess),
        stan_data$Y + 1))
  }
  
  fit <- mod$sample(
    data = stan_data,
    seed = seed_i, 
    chains = 6, 
    parallel_chains = 6,
    show_messages = TRUE, 
    show_exceptions = FALSE, 
    iter_warmup = 4000, 
    iter_sampling = 15000,
    save_warmup = FALSE,
    init = init_fun,
    thin=30
  )
  
  draws <- fit$draws()              # posterior draws (array)
  wb_rhat <- fit$summary()          # tibble with mean, sd, rhat, ess_*
  
  # draws <- as_draws_array(fit)
  # mcmc_trace(draws[, 5, , drop = FALSE], pars = c("X[4,1]"))
  # 
  # summarise_draws(draws[, -1, , drop = FALSE])
  # mcmc_trace(draws[, , , drop = FALSE], pars = c("X[4,120]"))
  
  # Return something small (paths, summary stats, etc.)
  file_name = paste0("stan_samples_nolag_noadm_wastewatersample_",waste_id,".rds")
 
 save(file = file_name, list = c("draws","wb_rhat"))
 
}

# record end time
end_time <- Sys.time()

# print runtime
print(end_time - start_time)

stopCluster(myCluster)


