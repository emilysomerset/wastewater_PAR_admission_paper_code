# parallel_foreach_thin_merge.R
library(data.table)
library(dplyr)
library(purrr)
library(foreach)
library(doParallel)

# --- extract_var_draws with thinning ---
extract_var_draws <- function(draws,
                              var,                   # e.g. "X", "beta", "alpha"
                              return = c("array", "long", "matrix"),
                              thin = 1) {
  return <- match.arg(return)
  
  if (!is.array(draws) || length(dim(draws)) != 3)
    stop("`draws` must be a 3D array with dims [iterations, chains, variables].")
  n_it  <- dim(draws)[1]
  n_ch  <- dim(draws)[2]
  vnames <- dimnames(draws)[[3]]
  if (is.null(vnames))
    stop("The 3rd dimension (variables) must have names: dimnames(draws)[[3]].")
  
  # thinning on iterations
  thin <- as.integer(thin)
  if (is.na(thin) || thin < 1) thin <- 1L
  it_idx <- seq(1, n_it, by = thin)
  if (length(it_idx) < n_it) {
    draws <- draws[it_idx, , , drop = FALSE]
    n_it <- dim(draws)[1]
  }
  
  pat_base <- paste0("^", var, "(?:\\[.*\\])?$")
  hits <- grepl(pat_base, vnames)
  if (!any(hits)) stop("No variable(s) matched `", var, "`.")
  sel_names <- vnames[hits]
  
  arr <- draws[, , hits, drop = FALSE]
  dim(arr) <- c(n_it * n_ch, dim(arr)[3])
  colnames(arr) <- sel_names
  n_draws_total <- nrow(arr)
  
  if (length(sel_names) == 1 && !grepl("\\[", sel_names)) {
    if (return == "matrix") return(arr)
    if (return == "long") {
      return(data.frame(
        variable = sel_names,
        draw = seq_len(n_draws_total),
        value = as.numeric(arr[, 1]),
        row.names = NULL
      ))
    }
    return(as.numeric(arr[, 1]))
  }
  
  if (return == "matrix") return(arr)
  
  if (return == "long") {
    out <- data.frame(
      variable = rep(sel_names, each = n_draws_total),
      draw = rep(seq_len(n_draws_total), times = length(sel_names)),
      value = as.vector(arr),
      row.names = NULL
    )
    return(out)
  }
  
  out <- arr
  dim(out) <- c(n_draws_total, length(sel_names))
  colnames(out) <- sel_names
  out
}



make_empty_acc <- function(n_grid) {
  acc <- vector("list", n_grid)
  for (i in seq_len(n_grid)) {
    acc[[i]] <- list(
      value = numeric(0),
      value_int = integer(0),
      Pi = numeric(0),
      Y = integer(0),
      age1 = numeric(0),
      age2 = numeric(0),
      age3 = numeric(0),
      age4 = numeric(0),
      X_int_3_4 = integer(0),
      X_cumsum = integer(0),
      X_cumsum_3_4 = integer(0),
      X_cumsum_delta = integer(0),
      X_cumsum_delta_3_4 = integer(0),
      total_draws = 0L
    )
  }
  acc
}

files <- list.files( pattern = "^stan_samples_nolag_noadm_wastewatersample_.*\\.rds$", full.names = TRUE)

J <- 185
I <- 4
grid <- expand.grid(age = 1:I, time = 1:J)
pairs <- paste0(grid$age, ",", grid$time)
n_grid <- nrow(grid)

thin <- 1  # thinning factor: 1 = no thinning
# Start cluster

myCluster <- makeCluster(20, # number of cores to use
                         type = "PSOCK")

registerDoParallel(myCluster)
print(getDoParWorkers())

# record start time
start_time <- Sys.time()

results <- foreach(f_idx = seq_len(length(files))) %dopar% {

  library(data.table)
  library(dplyr)
  library(purrr)
  library(foreach)
  library(doParallel)
  library(posterior)
  
  file <- files[f_idx]
  message(sprintf("Processing file %d / %d : %s", f_idx, length(files), basename(file)))
  
  load(file)
  conv_val <- max(wb_rhat$rhat)
  
  # if (conv_val > 1.05){
  #   j = 1
  #   while(conv_val >1.05 & j<= 6){
  #     conv_val = summarise_draws(draws[, -j, , drop = FALSE])%$% rhat %>% max()
  #     j = j + 1
  #   }
  #   if (conv_val <= 1.05){
  #     draws <- draws[,-(j-1), , drop = FALSE]
  #   }else{conv_val = max(wb_rhat$rhat)}
  # }
  
  if (conv_val < 1.05){  
    # prepare an empty accumulator for this file
    acc_file <- make_empty_acc(n_grid)
  
  # extract long formats with thinning
  X_long_all <- extract_var_draws(draws, var = "X", return = "long", thin = thin)
  Pi_long_all <- extract_var_draws(draws, var = "Pi", return = "long", thin = thin)
  C_long_all <- extract_var_draws(draws, var = "C", return = "long", thin = thin)

  # inside your foreach worker, after extracting X_long_all and Pi_long_all
  
  setDT(X_long_all)
  setDT(Pi_long_all)
  setDT(C_long_all)
  
  # extract time index from Pi; assume Pi[time]
  Pi_long_all[, time := as.integer(sub("^Pi\\[([0-9]+)\\]$", "\\1", variable))]
  pi_dt <- Pi_long_all[, .(draw, time, Pi = value)]
  setkey(pi_dt, draw, time)
  
  # extract time index from Pi; assume Pi[time]
  C_long_all[, c("age", "age_col") := tstrsplit(gsub("C\\[|\\]", "", variable), ",")]
  C_long_all[, age := as.integer(age)]
  C_long_all[, age_col := as.integer(age_col)]
  
  
  C_wide <- dcast(
    C_long_all,
    draw + age ~ age_col,
    value.var = "value"
  )
  
  setnames(C_wide, old = c("1","2","3","4"),
           new = c("age1","age2","age3","age4"))
  
  setkey(C_wide, draw)
  
  
  
  # extract age and time from X
  X_long_all[, c("age","time") := tstrsplit(sub("^X\\[(.*)\\]$", "\\1", variable), ",", type.convert = TRUE)]
  setkey(X_long_all, draw, time)
  
  # merge X and Pi by draw and time
  merged <- merge(X_long_all, pi_dt, by = c("draw", "time"))
  merged <- merge(merged, C_wide, by = c("draw", "age"))
  
  # sample Y ~ Binomial(X, Pi)
  # ensure X >= 0
  merged[, X_int := pmax(0, as.integer(ceiling(value)))]
  merged[, Y := rbinom(.N, size = X_int, prob = Pi)]
  merged[, X_cumsum_delta := cumsum(X_int)-X_int[1], by = .(draw,age)]
  merged[, X_cumsum := cumsum(X_int), by = .(draw,age)]
  
  #60 + calculations
  merged[age %in% 3:4, 
         X_int_3_4 := sum(X_int), 
         by = .(draw, time)]
  
  merged[age %in% 3:4, 
         X_cumsum_delta_3_4 := sum(X_cumsum_delta), 
         by = .(draw, time)]
  
  merged[age %in% 3:4, 
         X_cumsum_3_4 := sum(X_cumsum), 
         by = .(draw, time)]
  
  # now assign into accumulator
  for (i in seq_len(nrow(merged))) {
    age_i <- merged$age[i]
    time_j <- merged$time[i]
    key_grid <- match(paste0(age_i,",",time_j), pairs)
    
    acc_file[[key_grid]]$value <- c(acc_file[[key_grid]]$value, merged$value[i])
    acc_file[[key_grid]]$value_int <- c(acc_file[[key_grid]]$value_int, merged$X_int[i])
    acc_file[[key_grid]]$Pi <- c(acc_file[[key_grid]]$Pi, merged$Pi[i])
    acc_file[[key_grid]]$Y <- c(acc_file[[key_grid]]$Y, merged$Y[i])
    acc_file[[key_grid]]$age1 <- c(acc_file[[key_grid]]$age1, merged$age1[i])
    acc_file[[key_grid]]$age2 <- c(acc_file[[key_grid]]$age2, merged$age2[i])
    acc_file[[key_grid]]$age3 <- c(acc_file[[key_grid]]$age3, merged$age3[i])
    acc_file[[key_grid]]$age4 <- c(acc_file[[key_grid]]$age4, merged$age4[i])
    acc_file[[key_grid]]$X_int_3_4 <- c(acc_file[[key_grid]]$X_int_3_4, merged$X_int_3_4[i])
    acc_file[[key_grid]]$X_cumsum_delta <- c(acc_file[[key_grid]]$X_cumsum_delta, merged$X_cumsum_delta[i])
    acc_file[[key_grid]]$X_cumsum <- c(acc_file[[key_grid]]$X_cumsum, merged$X_cumsum[i])
    acc_file[[key_grid]]$X_cumsum_delta_3_4 <- c(acc_file[[key_grid]]$X_cumsum_delta_3_4, merged$X_cumsum_delta_3_4[i])
    acc_file[[key_grid]]$X_cumsum_3_4 <- c(acc_file[[key_grid]]$X_cumsum_3_4, merged$X_cumsum_3_4[i])
    acc_file[[key_grid]]$total_draws <- acc_file[[key_grid]]$total_draws + 1L
  }
  rm(draws, X_long_all, Pi_long_all)
  list(acc = acc_file, conv = conv_val)
  }
}

# stop cluster
stopCluster(myCluster)
message("Cluster stopped.")

# collect conv values
conv <- sapply(results, function(x) x$conv)
keep <- which(Reduce(c,lapply(conv, function(x)!is.null(x))))

results <- results[keep]

# merge per-file accumulators into final_acc
acc_list_per_file <- lapply(results, `[[`, "acc")
final_acc <- make_empty_acc(n_grid)

for (accf in acc_list_per_file) {
  for (i in seq_len(n_grid)) {
    if (length(accf[[i]]$value) == 0) next
    final_acc[[i]]$value <- c(final_acc[[i]]$value, accf[[i]]$value)
    final_acc[[i]]$value_int <- c(final_acc[[i]]$value_int, accf[[i]]$value_int)
    final_acc[[i]]$Y <- c(final_acc[[i]]$Y, accf[[i]]$Y)
    final_acc[[i]]$age1 <- c(final_acc[[i]]$age1, accf[[i]]$age1)
    final_acc[[i]]$age2 <- c(final_acc[[i]]$age2, accf[[i]]$age2)
    final_acc[[i]]$age3 <- c(final_acc[[i]]$age3, accf[[i]]$age3)
    final_acc[[i]]$age4 <- c(final_acc[[i]]$age4, accf[[i]]$age4)
    final_acc[[i]]$Pi <- c(final_acc[[i]]$Pi, accf[[i]]$Pi)
    final_acc[[i]]$X_int_3_4 <- c(final_acc[[i]]$X_int_3_4, accf[[i]]$X_int_3_4)
    final_acc[[i]]$X_cumsum_delta <- c(final_acc[[i]]$X_cumsum_delta, accf[[i]]$X_cumsum_delta)
    final_acc[[i]]$X_cumsum <- c(final_acc[[i]]$X_cumsum, accf[[i]]$X_cumsum)
    final_acc[[i]]$X_cumsum_delta_3_4 <- c(final_acc[[i]]$X_cumsum_delta_3_4, accf[[i]]$X_cumsum_delta_3_4)
    final_acc[[i]]$X_cumsum_3_4 <- c(final_acc[[i]]$X_cumsum_3_4, accf[[i]]$X_cumsum_3_4)
    final_acc[[i]]$total_draws <- final_acc[[i]]$total_draws + accf[[i]]$total_draws
  }
}

# --- compute final medians & quantiles ---
final_list <- vector("list", n_grid)
for (i in seq_len(n_grid)) {
  out <- list(age = grid$age[i], time = grid$time[i], total_draws = final_acc[[i]]$total_draws)
  
  nums <- c("value", "value_int", "Y","age1","age2","age3","age4", "Pi","X_int_3_4","X_cumsum_delta","X_cumsum_delta_3_4","X_cumsum","X_cumsum_3_4")
  for (nm in nums) {
    vec <- final_acc[[i]][[nm]]
    if (length(vec) == 0) {
      out[[paste0(nm, "_med")]] <- NA_real_
      out[[paste0(nm, "_q025")]] <- NA_real_
      out[[paste0(nm, "_q975")]] <- NA_real_
    } else {
      out[[paste0(nm, "_med")]] <- median(vec, na.rm = TRUE)
      out[[paste0(nm, "_q025")]] <- quantile(vec, probs = 0.025, na.rm = TRUE, names = FALSE)
      out[[paste0(nm, "_q05")]] <- quantile(vec, probs = 0.05, na.rm = TRUE, names = FALSE)
      out[[paste0(nm, "_q10")]] <- quantile(vec, probs = 0.1, na.rm = TRUE, names = FALSE)
      out[[paste0(nm, "_q90")]] <- quantile(vec, probs = 0.9, na.rm = TRUE, names = FALSE)
      out[[paste0(nm, "_q95")]] <- quantile(vec, probs = 0.95, na.rm = TRUE, names = FALSE)
      out[[paste0(nm, "_q975")]] <- quantile(vec, probs = 0.975, na.rm = TRUE, names = FALSE)
    }
  }
  final_list[[i]] <- as.data.frame(out, stringsAsFactors = FALSE)
}

final_summary <- bind_rows(final_list)
save(file = "results_nolag_toronto_thinned_noadm.RData", list = c("final_summary","conv"))

message("Done. Results saved to results_lagged_toronto.RData")
