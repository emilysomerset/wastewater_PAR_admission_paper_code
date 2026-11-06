
functions {
  real multinomial_log_approx(int Y, real X, real p) {
    if (X < Y) return negative_infinity();

    real x1 = Y;
    real x2 = X - Y;

    return lgamma(X + 1)
           - lgamma(x1 + 1)
           - lgamma(x2 + 1)
           + x1 * log(p)
           + x2 * log1m(p);
  }

  real poisson_log_approx(real X, real lambda) {
    if (lambda <= 0 || X <= -1)
      return negative_infinity();

    return X * log(lambda)
           - lambda
           - lgamma(X + 1);
  }
}

data {
  int<lower=1> I;
  int<lower=2> J;

  array[I, J] int<lower=0> Y;
  vector<lower=0>[J] ratio;

  vector<lower=0>[I] lower_X0;
  vector<lower=0>[I] mean_X0;
  vector<lower=0>[I] sd_X0;

  real<lower=0> phi_p;

}

parameters {
  // Global reporting trend
  real q1;
  real<lower=0> sigma_q;
  vector[J - 1] z_q;

  // Age-mixing matrix
  array[I] simplex[I] C;

  // Latent infections, including X_{i,0}
  vector[I] z_X0;
  vector[I * J] z_X_flat;
}

transformed parameters {
  vector[J] q;
  matrix<lower=0, upper=1>[I, J] Pi;

  vector<lower=0>[I] X_init;
  matrix<lower=0>[I, J] X;

  // -------------------------------------------
  // Global reporting probability random walk
  // -------------------------------------------

  q[1] = q1;

  for (j in 2:J) {
    q[j] = q[j - 1] + sigma_q * z_q[j - 1];
  }

  // Same reporting probability for every age group in a given week
  for (i in 1:I) {
    for (j in 1:J) {
      Pi[i, j] = inv_logit(q[j]);
    }
  }

  // -------------------------------------------
  // Initial latent infections X_{i,0}
  // -------------------------------------------

  for (i in 1:I) {
    real lower0 = lower_X0[i];
    real scale0 = sqrt(mean_X0[i] + 1);
    real a0 = (mean_X0[i] - lower0) / scale0;

    X_init[i] =
      lower0
      + scale0 * log1p_exp(a0 + z_X0[i]);
  }

  // -------------------------------------------
  // Latent infections X_{i,j}
  // -------------------------------------------

  for (j in 1:J) {
    for (i in 1:I) {
      int k = (i - 1) * J + j;

      real mixed_prev = 0;
      real lambda_ij;
      real lower_ij;
      real scale_ij;
      real a_ij;

      for (ip in 1:I) {
        if (j == 1) {
          mixed_prev += C[ip][i] * X_init[ip];
        } else {
          mixed_prev += C[ip][i] * X[ip, j - 1];
        }
      }

      lambda_ij = ratio[j] * mixed_prev;

      lower_ij = Y[i, j];
      scale_ij = sqrt(lambda_ij + 1);
      a_ij = (lambda_ij - lower_ij) / scale_ij;

      X[i, j] =
        lower_ij
        + scale_ij * log1p_exp(a_ij + z_X_flat[k]);
    }
  }
}

model {
  // -------------------------------------------
  // Reporting probability priors
  // -------------------------------------------

  q1 ~ normal(0, 1);
  z_q ~ std_normal();
  sigma_q ~ exponential(phi_p);

  // -------------------------------------------
  // Age-mixing matrix prior
  // -------------------------------------------

  for (i in 1:I) {
    C[i] ~ dirichlet(rep_vector(1.0, I));
  }

  // -------------------------------------------
  // Initial infections and Jacobians
  // -------------------------------------------

  for (i in 1:I) {
    real lower0 = lower_X0[i];
    real scale0 = sqrt(mean_X0[i] + 1);
    real a0 = (mean_X0[i] - lower0) / scale0;

    target += log(scale0)
              + log_inv_logit(a0 + z_X0[i]);

    target += normal_lpdf(
      X_init[i] | mean_X0[i], sd_X0[i]
    );
  }

  // -------------------------------------------
  // Infection dynamics and Jacobians
  // -------------------------------------------

  for (j in 1:J) {
    for (i in 1:I) {
      int k = (i - 1) * J + j;

      real mixed_prev = 0;
      real lambda_ij;
      real lower_ij;
      real scale_ij;
      real a_ij;

      for (ip in 1:I) {
        if (j == 1) {
          mixed_prev += C[ip][i] * X_init[ip];
        } else {
          mixed_prev += C[ip][i] * X[ip, j - 1];
        }
      }

      lambda_ij = ratio[j] * mixed_prev;

      lower_ij = Y[i, j];
      scale_ij = sqrt(lambda_ij + 1);
      a_ij = (lambda_ij - lower_ij) / scale_ij;

      // Jacobian adjustment
      target += log(scale_ij)
                + log_inv_logit(a_ij + z_X_flat[k]);

      // Infection process likelihood
      target += poisson_log_approx(X[i, j], lambda_ij);
    }
  }

  // -------------------------------------------
  // Observation likelihood
  // -------------------------------------------

  for (i in 1:I) {
    for (j in 1:J) {
      target += multinomial_log_approx(
        Y[i, j], X[i, j], Pi[i, j]
      );
    }
  }
}
