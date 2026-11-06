functions {
  real multinomial_log_approx(int A, int Y, real X, vector p) {
    if (X < Y || A < 0 || A > Y) return negative_infinity();
    real x1 = A;
    real x2 = Y - A;
    real x3 = X - Y;
    return lgamma(X + 1) - lgamma(x1 + 1) - lgamma(x2 + 1)
           - lgamma(x3 + 1) + x1 * log(p[1]) + x2 * log(p[2])
           + x3 * log(p[3]);
  }
  real poisson_log_approx(real X, real lambda) {
    if (lambda <= 0 || X <= -1) return negative_infinity();
    return X * log(lambda) - lambda - lgamma(X + 1);
  }
}

data {
  int<lower=1> I;
  int<lower=1> J;
  array[I, J] int<lower=0> Y;
  vector[J] ratio;
  vector<lower=0>[I] lower_X0;
  vector<lower=0>[I] mean_X0;
  vector<lower=0>[I] sd_X0;
  real<lower=0> phi_p;
  array[I, J] int<lower=0> A;
  real phi_tau;
  vector[I] tau_mean_prior;
}

parameters {
  // Unconstrained coordinates for the same infection variables.
  vector[I] z_X0;
  vector[I * J] z_X_flat;
  array[I] simplex[I] C;

  // Unchanged reporting RW and admissions AR(1) parameterization.
  vector[J] z_logit_pi;
  real<lower=0> sd_p;
  matrix[I, J] z_logit_tau;
  real<lower=0> sd_tau;
  vector[I] tau_mean;
  real<lower=0, upper=1> rho_tau;
}

transformed parameters {
  vector<lower=0>[I] X_init;
  matrix<lower=0>[I, J] X;
  vector[J] Pi;
  vector[J] logit_pi;
  matrix[I, J] tau;
  matrix[I, J] logit_tau;
  array[I, J] simplex[3] p;

  // Same RW as the original, including its initial distribution.
  logit_pi[1] = sd_p * z_logit_pi[1];
  for (j in 2:J)
    logit_pi[j] = logit_pi[j - 1] + sd_p * z_logit_pi[j];

  // Same AR(1) as the original, including its initial distribution.
  for (i in 1:I) {
    logit_tau[i, 1] = tau_mean[i] + sd_tau * z_logit_tau[i, 1];
    for (j in 2:J)
      logit_tau[i, j] = tau_mean[i]
        + (logit_tau[i, j - 1] - tau_mean[i]) * rho_tau
        + sd_tau * z_logit_tau[i, j];
  }

  for (j in 1:J) {
    Pi[j] = inv_logit(logit_pi[j]);
    for (i in 1:I) {
      tau[i, j] = inv_logit(logit_tau[i, j]);
      p[i, j][1] = tau[i, j];
      p[i, j][2] = (1 - tau[i, j]) * Pi[j];
      p[i, j][3] = (1 - tau[i, j]) * (1 - Pi[j]);
    }
  }

  // Softplus bijection: X_init[i] > lower_X0[i].
  for (i in 1:I) {
    real scale0 = sqrt(mean_X0[i] + 1);
    real a0 = (mean_X0[i] - lower_X0[i]) / scale0;
    X_init[i] = lower_X0[i]
      + scale0 * log1p_exp(a0 + z_X0[i]);
  }

  // Softplus bijection: X[i,j] > Y[i,j].
  // The original likelihood is zero when X[i,j] < Y[i,j].
  for (j in 1:J) {
    for (i in 1:I) {
      int k = (i - 1) * J + j;
      real mixed_prev = 0;
      real lambda_ij;
      real scale_ij;
      real a_ij;
      for (ip in 1:I) {
        if (j == 1)
          mixed_prev += C[ip][i] * X_init[ip];
        else
          mixed_prev += C[ip][i] * X[ip, j - 1];
      }
      lambda_ij = ratio[j] * mixed_prev;
      scale_ij = sqrt(lambda_ij + 1);
      a_ij = (lambda_ij - Y[i, j]) / scale_ij;
      X[i, j] = Y[i, j]
        + scale_ij * log1p_exp(a_ij + z_X_flat[k]);
    }
  }
}

model {
  // Exactly the original priors for the reporting and admissions processes.
  z_logit_pi ~ normal(0, 1);
  to_vector(z_logit_tau) ~ normal(0, 1);
  sd_p ~ exponential(phi_p);
  sd_tau ~ exponential(phi_tau);
  rho_tau ~ beta(1, 1);
  tau_mean ~ normal(tau_mean_prior, 0.5);

  for (i in 1:I)
    C[i] ~ dirichlet(rep_vector(1.0, I));

  // Original truncated-normal initial-infection prior, plus Jacobian.
  for (i in 1:I) {
    real scale0 = sqrt(mean_X0[i] + 1);
    real a0 = (mean_X0[i] - lower_X0[i]) / scale0;
    target += normal_lpdf(X_init[i] | mean_X0[i], sd_X0[i]);
    target += log(scale0) + log_inv_logit(a0 + z_X0[i]);
  }

  // Original infection likelihood, plus Jacobian for each X[i,j].
  for (j in 1:J) {
    for (i in 1:I) {
      int k = (i - 1) * J + j;
      real mixed_prev = 0;
      real lambda_ij;
      real scale_ij;
      real a_ij;
      for (ip in 1:I) {
        if (j == 1)
          mixed_prev += C[ip][i] * X_init[ip];
        else
          mixed_prev += C[ip][i] * X[ip, j - 1];
      }
      lambda_ij = ratio[j] * mixed_prev;
      scale_ij = sqrt(lambda_ij + 1);
      a_ij = (lambda_ij - Y[i, j]) / scale_ij;
      target += poisson_log_approx(X[i, j], lambda_ij);
      target += log(scale_ij) + log_inv_logit(a_ij + z_X_flat[k]);
    }
  }

  for (i in 1:I)
    for (j in 1:J)
      target += multinomial_log_approx(A[i, j], Y[i, j], X[i, j], p[i, j]);
}
