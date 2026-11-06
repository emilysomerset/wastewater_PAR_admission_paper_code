#!/usr/bin/env Rscript
# MCMC diagnostics for the admissions + cases model (common reporting RW).
# Each parameter's diagnostics are calculated within each of 500 conditional fits.
# Summaries pool parameter-fit combinations within each named group.
# The original draws and diagnostic objects are not modified.

library(dplyr)

root <- Sys.getenv(
  "ADMISSIONS_FIT_DIR",
  "/scratch/somerse5/wastewater_admissions_common_prior_sensitivity/fits"
)
output_dir <- Sys.getenv("DIAGNOSTICS_OUTPUT_DIR", dirname(root))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(root,
  pattern = "^admissions_common_revised_wastewatersample_[0-9]+\\.RData$",
  full.names = TRUE
))
if (length(files) != 500L) stop("Expected 500 fits; found ", length(files))

get_group <- function(variable) {
  case_when(
    variable == "lp__" ~ "Log posterior density",
    variable == "sd_p" ~ "Reporting RW SD",
    variable == "sd_tau" ~ "Hospitalization AR SD",
    variable == "rho_tau" ~ "Hospitalization AR correlation",
    grepl("^tau_mean\\[", variable) ~ "Baseline hospitalization logits",
    grepl("^(z_logit_tau|logit_tau|tau)\\[", variable) ~ "Hospitalization process",
    grepl("^(z_logit_pi|logit_pi|Pi)\\[", variable) ~ "Reporting process",
    grepl("^p\\[", variable) ~ "Multinomial observation probabilities",
    grepl("^C\\[", variable) ~ "Age-mixing weights",
    grepl("^(z_X0|X_init)\\[", variable) ~ "Initial infections",
    grepl("^(z_X_flat|X)\\[", variable) ~ "Latent infections",
    TRUE ~ "Other"
  )
}

summaries <- vector("list", length(files))
divergence_records <- vector("list", length(files))
for (i in seq_along(files)) {
  e <- new.env(parent = emptyenv())
  load(files[i], envir = e)
  if (!exists("wb_rhat", envir = e, inherits = FALSE))
    stop("Missing wb_rhat in ", files[i])
  x <- e$wb_rhat
  required <- c("variable", "rhat", "ess_bulk", "ess_tail")
  if (!all(required %in% names(x))) stop("Missing diagnostic columns in ", files[i])
  x <- x |>
    mutate(group = get_group(variable)) |>
    select(group, rhat, ess_bulk, ess_tail)
  if (any(x$group == "Other")) stop("Unclassified variables in ", files[i])
  if (any(!is.finite(as.matrix(x[, c("rhat", "ess_bulk", "ess_tail")]))))
    stop("Non-finite diagnostics in ", files[i])
  summaries[[i]] <- x

  waste_id <- as.integer(sub(".*wastewatersample_([0-9]+)\\.RData$", "\\1", basename(files[i])))
  div <- NA_real_
  if (exists("sampler_diagnostics", envir = e, inherits = FALSE)) {
    d <- e$sampler_diagnostics
    if (length(dim(d)) == 3L && "divergent__" %in% dimnames(d)[[3]])
      div <- sum(d[, , "divergent__"])
  }
  divergence_records[[i]] <- data.frame(waste_id = waste_id, divergences = div)
  rm(e)
  if (i %% 50L == 0L) message("Processed ", i, "/500 fits")
}

diagnostics <- bind_rows(summaries)
rhat_table <- diagnostics |>
  group_by(group) |>
  summarise(
    N = n(),
    median = median(rhat),
    q75 = as.numeric(quantile(rhat, .75)),
    max = max(rhat),
    prop_above_101 = mean(rhat > 1.01),
    prop_above_105 = mean(rhat > 1.05),
    .groups = "drop"
  )
ess_table <- diagnostics |>
  group_by(group) |>
  summarise(
    N = n(),
    bulk_min = min(ess_bulk),
    bulk_q25 = as.numeric(quantile(ess_bulk, .25)),
    bulk_median = median(ess_bulk),
    tail_min = min(ess_tail),
    tail_q25 = as.numeric(quantile(ess_tail, .25)),
    tail_median = median(ess_tail),
    .groups = "drop"
  )
divergences <- bind_rows(divergence_records)

write.csv(rhat_table, file.path(output_dir, "mcmc_diagnostics_admissions_rhat.csv"), row.names = FALSE)
write.csv(ess_table, file.path(output_dir, "mcmc_diagnostics_admissions_ess.csv"), row.names = FALSE)
write.csv(divergences, file.path(output_dir, "mcmc_diagnostics_admissions_divergences_by_fit.csv"), row.names = FALSE)
print(rhat_table, n = Inf)
print(ess_table, n = Inf)
message("Fits with Rhat > 1.01: ", sum(vapply(summaries, function(x) any(x$rhat > 1.01), logical(1))))
message("Total divergences: ", if (all(is.na(divergences$divergences))) "not available" else sum(divergences$divergences, na.rm = TRUE))
