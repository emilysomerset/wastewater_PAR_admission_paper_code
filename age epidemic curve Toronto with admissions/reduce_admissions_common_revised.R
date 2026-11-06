#!/usr/bin/env Rscript

root <- "/scratch/somerse5/wastewater_admissions_common_prior_sensitivity"

input_dir <- file.path(root, "fits")
output_dir <- file.path(root, "fits_reduced")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2917)
ids <- sample(1:3000, 500)

for (f in seq_along(ids)) {

  waste_id <- ids[f]

  filename <- sprintf(
    "admissions_common_revised_wastewatersample_%d.RData",
    waste_id
  )

  input <- file.path(input_dir, filename)
  output <- file.path(output_dir, filename)

  if (file.exists(output)) next
  if (!file.exists(input)) stop("Missing input: ", input)

  e <- new.env()
  load(input, envir = e)

  draws <- e$draws

  if (dim(draws)[2] != 4L || dim(draws)[1] != 2000L) {
    stop("Unexpected draw dimensions: ", filename)
  }

  # Select 250 draws independently within each chain
  set.seed(10000L + waste_id)

  idx <- sort(sample.int(2000L, 250L))

  draws <- draws[idx, , , drop = FALSE]

  wb_rhat <- e$wb_rhat
  sampler_diagnostics <- e$sampler_diagnostics

  tmp <- paste0(output, ".partial")

  save(
    draws, wb_rhat, sampler_diagnostics,
    file = tmp,
    compress = TRUE
  )

  if (!file.rename(tmp, output)) {
    stop("Could not save: ", output)
  }

  rm(e, draws, wb_rhat, sampler_diagnostics)
  gc(verbose = FALSE)

  if (f %% 50 == 0) {
    message("Reduced ", f, "/500 fits")
  }
}

message("Finished reducing all 500 fits.")
