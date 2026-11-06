#!/usr/bin/env Rscript
# MCMC diagnostics for case-only epidemic model, common reporting RW.
# Usage: Rscript mcmc_diagnostics_case_only.R [fit_directory] [output_directory]
# Diagnostics are computed per parameter within each conditional fit, then
# summarized over all parameter-fit combinations in each group.

args <- commandArgs(trailingOnly = TRUE)
fit_dir <- if (length(args) >= 1) args[[1]] else "/scratch/somerse5/wastewater_model/results_noadm_common_rw_long"
out_dir <- if (length(args) >= 2) args[[2]] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(fit_dir, pattern = "^stan_samples_common_rw_long_wastewatersample_[0-9]+\\.RData$", full.names = TRUE))
if (length(files) != 500L) stop("Expected 500 fits, found ", length(files))

group_of <- function(v) {
  g <- rep("Other", length(v))
  g[v == "lp__"] <- "Log posterior density"
  g[v == "sigma_q"] <- "Reporting RW SD"
  g[v == "q1"] <- "Initial reporting logit"
  g[grepl("^(z_q|q|Pi)\\[", v)] <- "Reporting process"
  g[grepl("^C\\[", v)] <- "Age-mixing weights"
  g[grepl("^(z_X0|X_init)\\[", v)] <- "Initial infections"
  g[grepl("^(z_X_flat|X)\\[", v)] <- "Latent infections"
  g
}

pieces <- vector("list", length(files))
for (i in seq_along(files)) {
  e <- new.env(parent = emptyenv())
  load(files[[i]], envir = e)
  x <- e$wb_rhat
  required <- c("variable", "rhat", "ess_bulk", "ess_tail")
  if (!all(required %in% names(x))) stop("Missing diagnostics in ", basename(files[[i]]))
  pieces[[i]] <- data.frame(group = group_of(x$variable), rhat = x$rhat,
                            ess_bulk = x$ess_bulk, ess_tail = x$ess_tail)
  rm(e, x)
  if (i %% 50L == 0L) message("Read ", i, "/500 fits")
}
d <- do.call(rbind, pieces)
if (any(d$group == "Other")) stop("Unclassified parameters detected; inspect Stan variable names before reporting")
if (any(!is.finite(d$rhat) | !is.finite(d$ess_bulk) | !is.finite(d$ess_tail))) stop("Nonfinite diagnostics")

spl <- split(d, d$group)
rhat <- do.call(rbind, lapply(names(spl), function(g) {
  x <- spl[[g]]$rhat
  data.frame(group = g, N = length(x), median = median(x), q75 = unname(quantile(x, .75)),
             max = max(x), prop_above_101 = mean(x > 1.01), prop_above_105 = mean(x > 1.05))
}))
ess <- do.call(rbind, lapply(names(spl), function(g) {
  x <- spl[[g]]$ess_bulk; y <- spl[[g]]$ess_tail
  data.frame(group = g, N = length(x), bulk_min = min(x), bulk_q25 = unname(quantile(x, .25)),
             bulk_median = median(x), tail_min = min(y), tail_q25 = unname(quantile(y, .25)),
             tail_median = median(y))
}))
rownames(rhat) <- NULL; rownames(ess) <- NULL
write.csv(rhat, file.path(out_dir, "mcmc_diagnostics_case_only_rhat.csv"), row.names = FALSE)
write.csv(ess, file.path(out_dir, "mcmc_diagnostics_case_only_ess.csv"), row.names = FALSE)
print(rhat, row.names = FALSE)
print(ess, row.names = FALSE)
