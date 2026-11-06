library(ggplot2)
library(dplyr)
library(tidyr)

root <- "/scratch/somerse5/wastewater_admissions_common_prior_sensitivity"
fit_dir <- file.path(root, "fits")

files <- list.files(
  fit_dir,
  pattern = "^admissions_common_revised_wastewatersample_[0-9]+\\.RData$",
  full.names = TRUE
)

stopifnot(length(files) == 500)

set.seed(2917)

parameters <- c(
  "rho_tau", "sd_tau", "sd_p",
  paste0("tau_mean[", 1:4, "]")
)

posterior_list <- vector("list", length(files))

for (f in seq_along(files)) {

  e <- new.env()
  load(files[f], envir = e)

  d <- e$draws
  vars <- dimnames(d)[[3]]

  if (!all(parameters %in% vars)) {
    stop("Missing parameters in ", basename(files[f]))
  }

  # Sample 1000 draws from the 8000 available
  n <- prod(dim(d)[1:2])
  idx <- sample.int(n, 1000)

  posterior_list[[f]] <- data.frame(
    rho_tau = as.vector(d[, , "rho_tau"])[idx],
    sd_tau = as.vector(d[, , "sd_tau"])[idx],
    sd_p = as.vector(d[, , "sd_p"])[idx],
    tau1 = plogis(as.vector(d[, , "tau_mean[1]"])[idx]),
    tau2 = plogis(as.vector(d[, , "tau_mean[2]"])[idx]),
    tau3 = plogis(as.vector(d[, , "tau_mean[3]"])[idx]),
    tau4 = plogis(as.vector(d[, , "tau_mean[4]"])[idx])
  )

  rm(e, d)

  if (f %% 50 == 0) {
    message("Processed ", f, "/500 fits")
  }
}

posterior <- bind_rows(posterior_list)

write.csv(
  posterior,
  file.path(root, "admissions_posterior_draws.csv"),
  row.names = FALSE
)



