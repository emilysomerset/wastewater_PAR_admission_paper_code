functions {
real multinomial_log_approx(int Y, real X, real p) {
      // validity check first
    if (X < Y) {
      return negative_infinity();
    }
    
  real x1 = Y;
  real x2 = X - Y;
  real log_pmf = lgamma(X + 1)
                 - lgamma(x1 + 1)
                 - lgamma(x2 + 1)
                 + x1 * log(p)
                 + x2 * log(1-p);
    return log_pmf;
}
real poisson_log_approx(real X, real lambda) {
    // Approximate log PMF of Poisson(Y | lambda)
    // using Stirling's approximation for log-factorial
    if (lambda <= 0||X <= -1) {
      return negative_infinity();
    }

    // Use lgamma(Y + 1) for log(Y!)
    real log_pmf = X * log(lambda) - lambda - lgamma(X + 1);
    return log_pmf;
  }
}

data {
  int<lower=1> I;                 // number of age groups
  int<lower=1> J;                 // number of weeks
  int<lower=0> J_future;          // Used for forecasting
  array[I, J] int<lower=0> Y;         // observed reported cases
  vector[J_future] ratio;                // R_j: external reproduction number
  vector<lower=0>[I] lower_X0;
  vector<lower=0>[I] mean_X0;     // prior mean for initial latent infections
  vector<lower=0>[I] sd_X0;       // prior SD for initial infections
  real<lower=0> phi_p;            // exponential rate for pi random walk
}

parameters {
  // Infection model
  vector<lower=lower_X0>[I] X_init;        // latent infections
  matrix<lower=0>[I, J] X;        // latent infections
  array[I] simplex[I] C;                // contact/transition weights (each row is a simplex)
  
  // Probability of mild reported 
  vector[J] z_logit_pi;      // standardized logit of shared reporting probability pi_j
  real <lower = 0> sd_p;
}

transformed parameters {
  // Infection model  
   // vector[I] nu= exp(log_nu);
  
  // Probability of mild reported 
   vector[J] Pi; 
   // real sd_p = exp(-0.5 * psi_p);
   vector[J] logit_pi;
   
  logit_pi[1] = sd_p * z_logit_pi[1];
  Pi[1] = inv_logit(logit_pi[1]);
  for (j in 2:J){
    logit_pi[j] =  logit_pi[j - 1] + sd_p * z_logit_pi[j];
    Pi[j] = inv_logit(logit_pi[j]);}
    
}

model {

  z_logit_pi ~ normal(0, 1);

  // Prior on initial latent infections
  for (i in 1:I)
    X_init[i] ~ normal(mean_X0[i], sd_X0[i]);
    
    // Infection dynamics
  for (j in 1:J) {
    for (i in 1:I) {
      real lambda_ij = 0;
      for (ip in 1:I){
        if (j ==1){
          lambda_ij += C[ip, i] * X_init[ip];
        }else{
        lambda_ij += C[ip, i] * X[ip, j - 1];}}
      lambda_ij *= ratio[j];
      target += poisson_log_approx(X[i, j],  lambda_ij);
    }
  }

for (i in 1:I) {
    for (j in 1:J) {
      target += multinomial_log_approx(Y[i, j], X[i, j], Pi[j]);
    }
  }
  
  
  // priors
  for (i in 1:I)
    C[i] ~ dirichlet(rep_vector(1.0, I));  // flat prior 
      
  sd_p ~ exponential(phi_p); 

}
