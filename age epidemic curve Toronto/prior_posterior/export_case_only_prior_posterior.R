root <- "/scratch/somerse5/wastewater_model/results_noadm_common_rw_long"

files <- list.files(
  root,
  pattern = "^stan_samples_common_rw_long_wastewatersample_[0-9]+\\.RData$",
  full.names = TRUE
)

stopifnot(length(files) == 500)

set.seed(2917)

posterior_list <- vector("list", length(files))

for (f in seq_along(files)) {

  e <- new.env()
  load(files[f], envir = e)

  d <- e$draws

  if (!"sigma_q" %in% dimnames(d)[[3]]) {
    stop("sigma_q not found in ", basename(files[f]))
  }

  # Select 1,000 of the 8,000 retained draws
  n <- prod(dim(d)[1:2])
  idx <- sample.int(n, 1000)

  posterior_list[[f]] <- data.frame(
    sigma_pi = as.vector(d[, , "sigma_q"])[idx]
  )

  rm(e, d)
  gc(verbose = FALSE)

  if (f %% 50 == 0) {
    message("Processed ", f, "/500 fits")
  }
}

posterior <- do.call(rbind, posterior_list)

write.csv(
  posterior,
  file.path(root, "case_only_posterior_draws.csv"),
  row.names = FALSE
)

cat("Posterior draws:", nrow(posterior), "\n")
cat("Posterior median:", median(posterior$sigma_pi), "\n")
cat("Posterior 95% interval:",
    quantile(posterior$sigma_pi, c(0.025, 0.975)), "\n")
