
#!/usr/bin/env Rscript

# Pool posterior draws from the 500 conditional admissions + cases model fits.
#
# Each fit corresponds to one wastewater posterior sample and contains
# 1,000 retained draws (250 from each of four MCMC chains).
#
# The script calculates pooled posterior quantiles for latent infections,
# reporting probabilities, hospitalization probabilities, posterior
# predictive observations, and cumulative infections.
#
# Weeks are processed in batches to limit memory usage.
# Cumulative infections are calculated within each posterior draw
# before posterior quantiles are obtained.


# ------------------------------------------------------------------
# Directories and settings
# ------------------------------------------------------------------

root <- "/scratch/somerse5/wastewater_admissions_common_prior_sensitivity"
fit_dir <- file.path(root, "fits_reduced")

output <- file.path(root, "admissions_common_revised_pooled_reduced.RData")

# Reproduce the 500 wastewater samples used in model fitting
set.seed(2917)
ids <- sample(1:3000, 500)

weeks_per_batch <- 2L

files <- file.path(
  fit_dir,
  sprintf("admissions_common_revised_wastewatersample_%d.RData", ids)
)

if (any(!file.exists(files))) {
  stop(
    "Missing fit files: ",
    paste(head(basename(files[!file.exists(files)]), 5), collapse = ", ")
  )
}

# Posterior quantiles to calculate
probs <- c(.025, .05, .10, .50, .90, .95, .975)
suffix <- c("q025", "q05", "q10", "med", "q90", "q95", "q975")


# ------------------------------------------------------------------
# Determine model dimensions and number of posterior draws
# ------------------------------------------------------------------

e <- new.env(parent = emptyenv())
load(files[1], envir = e)

vars0 <- dimnames(e$draws)[[3]]

I <- length(grep("^tau_mean\\[[0-9]+\\]$", vars0))
J <- length(grep("^Pi\\[[0-9]+\\]$", vars0))

n_per <- prod(dim(e$draws)[1:2])

rm(e)

if (I != 4L || J < 1L) {
  stop("Unexpected model dimensions")
}

n_total <- n_per * length(files)

message(
  "Pooling ", length(files), " fits x ", n_per,
  " draws = ", n_total,
  "; I=", I, " J=", J
)


# ------------------------------------------------------------------
# Helper function for extracting posterior draws
# ------------------------------------------------------------------

# Extract specified Stan variables from a draws array with dimensions
# iteration x chain x variable.
#
# Returns a matrix with one row per posterior draw and one column
# per requested variable, preserving alignment between parameters.

extract <- function(draws, names) {
  
  existing <- dimnames(draws)[[3]]
  
  if (!all(names %in% existing)) {
    stop(
      "Missing parameter: ",
      paste(head(setdiff(names, existing), 4), collapse = ", ")
    )
  }
  
  a <- draws[, , names, drop = FALSE]
  
  matrix(
    as.numeric(a),
    nrow = prod(dim(a)[1:2]),
    ncol = length(names)
  )
}


# ------------------------------------------------------------------
# Variables to summarize
# ------------------------------------------------------------------

# Combined 60+ summaries are recorded for ages 3 and 4.
# Values for ages 1 and 2 are intentionally left as NA.

sixty <- c(
  "X_int_3_4",
  "X_cumsum_3_4",
  "X_cumsum_delta_3_4"
)

vars <- c(
  "value",
  "value_int",
  "Pi",
  "tau",
  "pY",
  "Y",
  "A",
  "tau_mean",
  "tau_mean_logit",
  "rho_tau",
  paste0("age", 1:I),
  "X_int_3_4",
  "X_cumsum",
  "X_cumsum_3_4",
  "X_cumsum_delta",
  "X_cumsum_delta_3_4"
)

# Note: tau_mean is the baseline hospitalization logit.
# The existing output variable tau_mean_logit contains
# plogis(tau_mean), i.e., the baseline hospitalization probability.
# Its original name is retained for compatibility with plotting code.


# ------------------------------------------------------------------
# Convergence diagnostics and scalar parameters
# ------------------------------------------------------------------

# Diagnostics were calculated from the original full MCMC chains
# and retained in the reduced posterior files.

conv <- data.frame(
  waste_id = ids,
  max_rhat = NA_real_,
  min_bulk_ess = NA_real_,
  min_tail_ess = NA_real_,
  divergences = NA_integer_
)

# Scalar model parameters to summarize across all wastewater fits
scalar <- c(
  "sd_p",
  "sd_tau",
  "rho_tau",
  sprintf("tau_mean[%d]", seq_len(I))
)

scalar_pool <- setNames(
  lapply(scalar, function(z) numeric(n_total)),
  scalar
)

# Store posterior summaries for each batch of weeks
all_summary <- vector("list", ceiling(J / weeks_per_batch))


# ------------------------------------------------------------------
# Pool posterior draws in batches of weeks
# ------------------------------------------------------------------

for (batch in seq_along(all_summary)) {
  
  weeks <- seq.int(
    (batch - 1L) * weeks_per_batch + 1L,
    min(batch * weeks_per_batch, J)
  )
  
  nw <- length(weeks)
  nc <- I * nw
  
  message(
    "Batch ", batch, "/", length(all_summary),
    ": weeks ", min(weeks), "-", max(weeks)
  )
  
  # Allocate matrices for the current batch.
  # Rows correspond to pooled posterior draws.
  # Columns correspond to age-week combinations.
  
  pool <- setNames(
    lapply(
      vars,
      function(z) matrix(NA_real_, nrow = n_total, ncol = nc)
    ),
    vars
  )
  
  # --------------------------------------------------------------
  # Process each conditional posterior fit
  # --------------------------------------------------------------
  
  for (f in seq_along(files)) {
    
    e <- new.env(parent = emptyenv())
    load(files[f], envir = e)
    
    d <- e$draws
    
    if (prod(dim(d)[1:2]) != n_per) {
      stop("Unexpected draw count in ", files[f])
    }
    
    # Rows occupied by this wastewater fit in the pooled matrices
    rows <- (f - 1L) * n_per + seq_len(n_per)
    
    # Extract diagnostics and scalar parameters only once
    if (batch == 1L) {
      
      sm <- e$wb_rhat
      
      conv$max_rhat[f] <- max(sm$rhat, na.rm = TRUE)
      conv$min_bulk_ess[f] <- min(sm$ess_bulk, na.rm = TRUE)
      conv$min_tail_ess[f] <- min(sm$ess_tail, na.rm = TRUE)
      
      if (
        exists("sampler_diagnostics", envir = e) &&
        "divergent__" %in% dimnames(e$sampler_diagnostics)[[3]]
      ) {
        conv$divergences[f] <- sum(
          e$sampler_diagnostics[, , "divergent__"]
        )
      }
      
      for (nm in scalar) {
        scalar_pool[[nm]][rows] <- extract(d, nm)[, 1L]
      }
    }
    
    # Common reporting probability and AR(1) correlation
    Pi <- extract(d, sprintf("Pi[%d]", weeks))
    rho <- extract(d, "rho_tau")[, 1L]
    
    age_values <- vector("list", I)
    
    # ------------------------------------------------------------
    # Calculate age-specific quantities
    # ------------------------------------------------------------
    
    for (age in seq_len(I)) {
      
      # Continuous latent infection counts
      X <- extract(
        d,
        sprintf("X[%d,%d]", age, seq_len(J))
      )
      
      # Convert latent infections to nonnegative integer counts
      Xi <- pmax(ceiling(X), 0)
      
      # Cumulative infections are calculated within each draw
      cum <- Xi
      
      if (J > 1L) {
        for (j in 2:J) {
          cum[, j] <- cum[, j - 1L] + Xi[, j]
        }
      }
      
      # Age-specific hospitalization probabilities
      tau <- extract(
        d,
        sprintf("tau[%d,%d]", age, weeks)
      )
      
      # ----------------------------------------------------------
      # Posterior predictive simulation
      # ----------------------------------------------------------
      
      # Three mutually exclusive observation categories:
      # 1. Hospital admission
      # 2. Reported non-hospitalized infection
      # 3. Unreported non-hospitalized infection
      #
      # Admissions are generated first, followed by reported
      # non-hospitalized infections among the remaining infections.
      #
      # A = simulated admissions
      # Y = simulated admissions + reported non-admissions
      
      set.seed(100000L + f + 1000L * batch + 100000L * age)
      
      aa <- matrix(
        rbinom(
          length(tau),
          size = as.vector(Xi[, weeks, drop = FALSE]),
          prob = as.vector(tau)
        ),
        nrow = n_per
      )
      
      yy_non <- matrix(
        rbinom(
          length(tau),
          size = as.vector(Xi[, weeks, drop = FALSE] - aa),
          prob = as.vector(Pi)
        ),
        nrow = n_per
      )
      
      # Baseline hospitalization logit and age-mixing weights
      tm <- extract(
        d,
        sprintf("tau_mean[%d]", age)
      )[, 1L]
      
      C <- extract(
        d,
        sprintf("C[%d,%d]", age, seq_len(I))
      )
      
      # Store quantities for this age group
      age_values[[age]] <- list(
        value = X[, weeks, drop = FALSE],
        value_int = Xi[, weeks, drop = FALSE],
        Pi = Pi,
        tau = tau,
        pY = (1 - tau) * Pi,
        Y = aa + yy_non,
        A = aa,
        tau_mean = matrix(tm, nrow = n_per, ncol = nw),
        tau_mean_logit = matrix(
          plogis(tm), nrow = n_per, ncol = nw
        ),
        rho_tau = matrix(rho, nrow = n_per, ncol = nw),
        X_cumsum = cum[, weeks, drop = FALSE],
        X_cumsum_delta = sweep(
          cum[, weeks, drop = FALSE],
          1L,
          Xi[, 1L],
          "-"
        ),
        C = C
      )
    }
    
    # ------------------------------------------------------------
    # Add age-specific draws to the pooled matrices
    # ------------------------------------------------------------
    
    for (age in seq_len(I)) {
      
      cols <- (age - 1L) * nw + seq_len(nw)
      v <- age_values[[age]]
      
      for (nm in c(
        "value", "value_int", "Pi", "tau", "pY",
        "Y", "A", "tau_mean", "tau_mean_logit",
        "rho_tau", "X_cumsum", "X_cumsum_delta"
      )) {
        pool[[nm]][rows, cols] <- v[[nm]]
      }
      
      # Age-mixing weights
      for (k in seq_len(I)) {
        pool[[paste0("age", k)]][rows, cols] <- v$C[, k]
      }
      
      # Combined quantities for ages 60+
      if (age >= 3L) {
        
        pool$X_int_3_4[rows, cols] <-
          age_values[[3]]$value_int +
          age_values[[4]]$value_int
        
        pool$X_cumsum_3_4[rows, cols] <-
          age_values[[3]]$X_cumsum +
          age_values[[4]]$X_cumsum
        
        pool$X_cumsum_delta_3_4[rows, cols] <-
          age_values[[3]]$X_cumsum_delta +
          age_values[[4]]$X_cumsum_delta
      }
    }
    
    rm(e, d, age_values)
    
    if (f %% 50L == 0L) {
      message("  Read ", f, "/", length(files), " fits")
    }
  }
  
  # --------------------------------------------------------------
  # Calculate posterior quantiles for this batch
  # --------------------------------------------------------------
  
  # Column order: all weeks for age 1, then age 2, etc.
  grid <- expand.grid(time = weeks, age = seq_len(I))
  grid <- grid[, c("age", "time")]
  
  grid$total_draws <- n_total
  
  for (nm in vars) {
    
    q <- vapply(
      seq_len(nc),
      function(k) {
        
        if (nm %in% sixty && grid$age[k] < 3L) {
          return(rep(NA_real_, length(probs)))
        }
        
        quantile(
          pool[[nm]][, k],
          probs = probs,
          names = FALSE
        )
      },
      numeric(length(probs))
    )
    
    for (k in seq_along(probs)) {
      grid[[paste0(nm, "_", suffix[k])]] <- q[k, ]
    }
    
    pool[[nm]] <- NULL
  }
  
  all_summary[[batch]] <- grid
  
  rm(pool, grid)
  gc(verbose = FALSE)
  
  # Save a checkpoint after each completed batch
  checkpoint_dir <- file.path(root, "checkpoints")
  
  dir.create(
    checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    all_summary[[batch]],
    file.path(
      checkpoint_dir,
      sprintf(
        "weeks_%03d_%03d.rds",
        min(weeks),
        max(weeks)
      )
    )
  )
}


# ------------------------------------------------------------------
# Summarize scalar model parameters
# ------------------------------------------------------------------

parameter_summary <- do.call(
  rbind,
  lapply(scalar, function(nm) {
    
    q <- quantile(
      scalar_pool[[nm]],
      probs = probs,
      names = FALSE
    )
    
    z <- data.frame(
      variable = nm,
      total_draws = n_total
    )
    
    for (k in seq_along(probs)) {
      z[[suffix[k]]] <- q[k]
    }
    
    z
  })
)


# ------------------------------------------------------------------
# Combine age-week summaries
# ------------------------------------------------------------------

final_summary <- do.call(rbind, all_summary)

final_summary <- final_summary[
  order(final_summary$time, final_summary$age),
]

row.names(final_summary) <- NULL


# ------------------------------------------------------------------
# Save pooled results
# ------------------------------------------------------------------

# Save to a temporary file first so that an interrupted write
# does not overwrite an existing completed result.

tmp <- paste0(output, ".partial")

save(
  final_summary,
  parameter_summary,
  conv,
  file = tmp
)

if (!file.rename(tmp, output)) {
  stop("Could not publish result")
}

message(
  "Saved ", output,
  "; ", nrow(final_summary), " rows"
)
