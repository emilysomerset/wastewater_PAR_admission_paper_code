compute_post_fun <- function (samps, global_samps = NULL, knots, refined_x, p, degree = 0) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps = matrix(0, nrow = p, ncol = ncol(samps))
  }
  if (nrow(global_samps) != p | nrow(samps) != (length(knots) - 
                                                1)) {
    return(message("Error: Incorrect dimension of global_samps or samps. Check whether the choice of p or the choice of knots are consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  X = global_poly_helper(refined_x, p = p)
  X <- as.matrix(X[, 1:(p - degree)])
  for (i in 1:ncol(X)) {
    X[, i] <- (factorial(i + degree - 1)/factorial(i - 1)) * 
      X[, i]
  }
  B = as(local_poly_helper(knots, refined_x = refined_x, p = (p - 
                                                                degree)), "dgTMatrix")
  fitted_samps_deriv <- X %*% global_samps[(1 + degree):p, 
  ] + B %*% samps
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}

process_results <- function(df_full, tmbdat, samps1, polyOrder){
  
  P = as.matrix(tmbdat$P)
  X = as.matrix(tmbdat$X)
  knots = tmbdat$knots
  
  coefsamps1 <- samps1$samps[1:ncol(P),]
  global_samps1 <- samps1$samps[(ncol(P) + 1):(ncol(P) + ncol(X)),]
  
  v <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                         knots = knots, 
                         refined_x = df_full$t,
                         p = polyOrder, degree = 0)
  
  vderiv <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                              knots = knots, 
                              refined_x = df_full$t,
                              p = polyOrder, degree = 1)
  
  ## Log Ospline 
  df_full$v <- as.numeric(apply(v[,-1], MARGIN=1,median))
  df_full$v_upr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.975))
  df_full$v_lwr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.025))
  
  df_full$v_deriv <- as.numeric(apply(vderiv[,-1], MARGIN=1,median))
  df_full$v_deriv_upr <- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.975))
  df_full$v_deriv_lwr<- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.025))
  
  ## Ospline
  df_full$exp_v <- as.numeric(apply(exp(v[,-1]), MARGIN=1,median))
  df_full$exp_v_upr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.975))
  df_full$exp_v_lwr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.025))
  
  df_full$exp_v_deriv <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,median))
  df_full$exp_v_deriv_upr <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.975))
  df_full$exp_v_deriv_lwr<- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.025))
  
  # Reproduction. number
  df_full$inst_repro = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,median))
  df_full$inst_repro_lwr = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,quantile, 0.025))
  df_full$inst_repro_upr = as.numeric(apply(exp(vderiv[,-1]), MARGIN=1,quantile, 0.975))

  return(list(df_full=df_full))
}