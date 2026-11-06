
#!/usr/bin/env Rscript

# Pool posterior draws from the 500 conditional case-only model fits.
#
# Each fit corresponds to one wastewater posterior sample and contains
# 1,000 preselected draws (250 from each of four MCMC chains).
#
# The script calculates pooled posterior quantiles for latent infections,
# reporting probabilities, posterior predictive case counts, age-mixing
# weights, and cumulative infections.
#
# Weeks are processed in batches to limit memory usage.
# Cumulative infections are calculated within each posterior draw
# before posterior quantiles are obtained.


# ------------------------------------------------------------------
# Directories and settings
# ------------------------------------------------------------------

result_dir <- "/scratch/somerse5/wastewater_model/results_noadm_common_rw_long_reduced"

output_file <- "/scratch/somerse5/wastewater_model/results_nolag_toronto_pooled_common_rw.RData"

# Model dimensions: four age groups and 185 weeks
I <- 4L
J <- 185L

# Number of weeks processed simultaneously
weeks_per_batch <- 10L

# Posterior quantiles
probs <- c(.025, .05, .10, .50, .90, .95, .975)
prob_suffix <- c("q025", "q05", "q10", "med", "q90", "q95", "q975")

# Locate the 500 reduced posterior files
files <- sort(list.files(
  result_dir,
  pattern = "^stan_samples_common_rw_wastewatersample_[0-9]+\\.RData$",
  full.names = TRUE
))

if (length(files) != 500L) {
  stop("Expected 500 result files; found ", length(files))
}


# ------------------------------------------------------------------
# Determine number of posterior draws
# ------------------------------------------------------------------

first <- new.env(parent = emptyenv())
load(files[1], envir = first)

n_per_fit <- prod(dim(first$draws)[1:2])

rm(first)

if (n_per_fit != 1000L) {
  stop("Expected 1,000 draws per fit; found ", n_per_fit)
}

n_total <- length(files) * n_per_fit

message(
  "Pooling ", length(files), " fits x ", n_per_fit,
  " draws = ", n_total, " draws per age-week."
)


# ------------------------------------------------------------------
# Helper function for extracting posterior draws
# ------------------------------------------------------------------

# Extract specified Stan variables from an array with dimensions
# iteration x chain x variable.
#
# Returns a matrix with one row per posterior draw and one column
# per requested variable, preserving alignment between parameters.

extract_matrix <- function(draws, names) {
  
  existing <- dimnames(draws)[[3]]
  missing <- setdiff(names, existing)
  
  if (length(missing)) {
    stop(
      "Missing Stan variable(s): ",
      paste(head(missing, 5), collapse = ", ")
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

vars <- c(
  "value",
  "value_int",
  "Pi",
  "q",
  "Y",
  paste0("age", 1:I),
  "X_int_3_4",
  "X_cumsum",
  "X_cumsum_3_4",
  "X_cumsum_delta",
  "X_cumsum_delta_3_4"
)

# Combined 60+ summaries are recorded for ages 3 and 4.
# Values for ages 1 and 2 are intentionally left as NA.

sixty_vars <- c(
  "X_int_3_4",
  "X_cumsum_3_4",
  "X_cumsum_delta_3_4"
)


# ------------------------------------------------------------------
# Convergence diagnostics and scalar parameters
# ------------------------------------------------------------------

# Diagnostics were calculated from the original full MCMC chains
# and retained in the reduced posterior files.

conv <- data.frame(
  waste_id = as.integer(sub(
    ".*wastewatersample_([0-9]+)\\.RData$",
    "\\1",
    basename(files)
  )),
  max_rhat = NA_real_,
  min_bulk_ess = NA_real_,
  min_tail_ess = NA_real_
)

# Scalar model parameters summarized across all wastewater fits
scalar_names <- c("q1", "sigma_q")

scalar_pool <- setNames(
  lapply(scalar_names, function(x) numeric(n_total)),
  scalar_names
)

# Store posterior summaries for each batch
all_summaries <- vector(
  "list",
  ceiling(J / weeks_per_batch)
)


# ------------------------------------------------------------------
# Pool posterior draws in batches of weeks
# ------------------------------------------------------------------

for (batch in seq_along(all_summaries)) {
  
  weeks <- seq.int(
    (batch - 1L) * weeks_per_batch + 1L,
    min(batch * weeks_per_batch, J)
  )
  
  n_weeks <- length(weeks)
  n_cells <- I * n_weeks
  
  message(
    "Batch ", batch, "/", length(all_summaries),
    ": weeks ", min(weeks), "-", max(weeks),
    " (", Sys.time(), ")"
  )
  
  # Allocate matrices for the current batch.
  # Rows correspond to pooled posterior draws.
  # Columns correspond to age-week combinations.
  
  pool <- setNames(
    lapply(
      vars,
      function(z) matrix(
        NA_real_,
        nrow = n_total,
        ncol = n_cells
      )
    ),
    vars
  )
  
  # --------------------------------------------------------------
  # Process each conditional posterior fit
  # --------------------------------------------------------------
  
  for (f in seq_along(files)) {
    
    e <- new.env(parent = emptyenv())
    load(files[f], envir = e)
    
    draws <- e$draws
    
    if (prod(dim(draws)[1:2]) != n_per_fit) {
      stop("Unexpected draw count in ", files[f])
    }
    
    # Rows occupied by this wastewater fit in the pooled matrices
    rows <- (f - 1L) * n_per_fit + seq_len(n_per_fit)
    
    # Extract convergence diagnostics and scalar parameters once
    if (batch == 1L) {
      
      s <- e$wb_rhat
      
      conv$max_rhat[f] <- max(s$rhat, na.rm = TRUE)
      conv$min_bulk_ess[f] <- min(s$ess_bulk, na.rm = TRUE)
      conv$min_tail_ess[f] <- min(s$ess_tail, na.rm = TRUE)
      
      for (nm in scalar_names) {
        scalar_pool[[nm]][rows] <-
          extract_matrix(draws, nm)[, 1L]
      }
    }
    
    # Reporting logits for the current weeks
    q_values <- extract_matrix(
      draws,
      sprintf("q[%d]", weeks)
    )
    
    age_values <- vector("list", I)
    
    # ------------------------------------------------------------
    # Calculate age-specific quantities
    # ------------------------------------------------------------
    
    for (age in seq_len(I)) {
      
      # Continuous latent infection counts
      X <- extract_matrix(
        draws,
        sprintf("X[%d,%d]", age, seq_len(J))
      )
      
      # Convert to nonnegative integer infection counts
      Xi <- pmax(ceiling(X), 0)
      
      # Calculate cumulative infections within each posterior draw
      cumulative <- Xi
      
      if (J > 1L) {
        for (j in 2:J) {
          cumulative[, j] <-
            cumulative[, j - 1L] + Xi[, j]
        }
      }
      
      # Cumulative infections excluding the initial week
      initial <- Xi[, 1L]
      
      delta <- sweep(
        cumulative[, weeks, drop = FALSE],
        1L,
        initial,
        "-"
      )
      
      # Age-specific reporting probabilities
      P <- extract_matrix(
        draws,
        sprintf("Pi[%d,%d]", age, weeks)
      )
      
      if (any(!is.finite(P)) || any(P < 0 | P > 1)) {
        stop("Invalid reporting probability in ", files[f])
      }
      
      # ----------------------------------------------------------
      # Posterior predictive reported case counts
      # ----------------------------------------------------------
      
      # Simulate reported cases conditional on latent infections
      # and the reporting probability for each posterior draw.
      #
      # Use a reproducible seed for each fit, batch, and age group.
      
      set.seed(
        100000L + f + 1000L * batch + 100000L * age
      )
      
      Y <- matrix(
        rbinom(
          length(P),
          size = as.vector(Xi[, weeks, drop = FALSE]),
          prob = as.vector(P)
        ),
        nrow = n_per_fit,
        ncol = n_weeks
      )
      
      # Age-mixing weights
      C <- extract_matrix(
        draws,
        sprintf("C[%d,%d]", age, seq_len(I))
      )
      
      age_values[[age]] <- list(
        value = X[, weeks, drop = FALSE],
        value_int = Xi[, weeks, drop = FALSE],
        Pi = P,
        q = q_values,
        Y = Y,
        X_cumsum = cumulative[, weeks, drop = FALSE],
        X_cumsum_delta = delta,
        C = C
      )
    }
    
    # ------------------------------------------------------------
    # Add age-specific draws to pooled matrices
    # ------------------------------------------------------------
    
    for (age in seq_len(I)) {
      
      cols <- (age - 1L) * n_weeks + seq_len(n_weeks)
      v <- age_values[[age]]
      
      for (nm in c(
        "value", "value_int", "Pi", "q",
        "Y", "X_cumsum", "X_cumsum_delta"
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
    
    rm(e, draws, age_values)
    
    if (f %% 50L == 0L) {
      message("  Read ", f, "/", length(files), " files")
    }
  }
  
  # --------------------------------------------------------------
  # Calculate posterior quantiles for the current batch
  # --------------------------------------------------------------
  
  # Column order: all weeks for age 1, then age 2, etc.
  grid <- expand.grid(time = weeks, age = seq_len(I))
  grid <- grid[, c("age", "time")]
  
  grid$total_draws <- n_total
  
  for (nm in vars) {
    
    q <- vapply(
      seq_len(n_cells),
      function(k) {
        
        if (nm %in% sixty_vars && grid$age[k] < 3L) {
          return(rep(NA_real_, length(probs)))
        }
        
        quantile(
          pool[[nm]][, k],
          probs = probs,
          names = FALSE,
          na.rm = FALSE
        )
      },
      numeric(length(probs))
    )
    
    for (p in seq_along(probs)) {
      grid[[paste0(nm, "_", prob_suffix[p])]] <- q[p, ]
    }
    
    pool[[nm]] <- NULL
    gc(verbose = FALSE)
  }
  
  all_summaries[[batch]] <- grid
  
  rm(pool, grid)
  gc(verbose = FALSE)
  
  # Save a checkpoint after each completed batch
  checkpoint_dir <- file.path(
    result_dir,
    "pooled_checkpoints"
  )
  
  dir.create(
    checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    all_summaries[[batch]],
    file.path(
      checkpoint_dir,
      sprintf(
        "pooled_checkpoint_weeks_%03d_%03d.rds",
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
  lapply(scalar_names, function(nm) {
    
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
      z[[prob_suffix[k]]] <- q[k]
    }
    
    z
  })
)


# ------------------------------------------------------------------
# Combine age-week summaries
# ------------------------------------------------------------------

final_summary <- do.call(rbind, all_summaries)

final_summary <- final_summary[
  order(final_summary$time, final_summary$age),
]

row.names(final_summary) <- NULL


# ------------------------------------------------------------------
# Save pooled results
# ------------------------------------------------------------------

# Write to a temporary file before replacing the final output.
# This prevents an interrupted save from overwriting a valid result.

tmp_output <- paste0(output_file, ".partial")

save(
  final_summary,
  parameter_summary,
  conv,
  file = tmp_output
)

if (!file.rename(tmp_output, output_file)) {
  stop("Could not move completed output to ", output_file)
}

message(
  "Saved final_summary (", nrow(final_summary),
  " rows), parameter_summary and conv to ", output_file
)
