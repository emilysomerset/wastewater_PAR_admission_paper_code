functions {
real multinomial_log_approx(int A, int Y, real X, vector p) {
      // validity check first
    if (X < Y || A < 0 || A > Y) {
      return negative_infinity();
    }
    
  real x1 = A;
  real x2 = Y - A;
  real x3 = X - Y;
  real log_pmf = lgamma(X + 1)
                 - lgamma(x1 + 1)
                 - lgamma(x2 + 1)
                 - lgamma(x3 + 1)
                 + x1 * log(p[1])
                 + x2 * log(p[2])
                 + x3 * log(p[3]);
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
  
  // Admissions
  array[I, J] int<lower=0> A;            // observed hospital admissions for COVID-19
  real phi_tau;                   // exponential rate for tau autoregressive
    vector[I] tau_mean_prior;
}

parameters {
  // Infection model
  vector<lower=lower_X0>[I] X_init;        // latent infections
  matrix<lower=0>[I, J] X;        // latent infections
  array[I] simplex[I] C;                // contact/transition weights (each row is a simplex)
  
  // Probability of mild reported 
  vector[J] z_logit_pi;      // standardized logit of shared reporting probability pi_j
  // real<lower=0> psi_p;            // log sd of pi RW(1)
  real <lower = 0> sd_p;
  // real<lower=0, upper=1> rho_pi;       // Autoregressive correlation parameter
  
  // Probability of severe reported
  matrix[I,J] z_logit_tau; 
  // real<lower=0> psi_tau;             // AR(1) scale
  real <lower = 0> sd_tau;
  vector[I] tau_mean;
  real<lower=0, upper=1> rho_tau;       // Autoregressive correlation parameter
}

transformed parameters {
  // Infection model  
   // vector[I] nu= exp(log_nu);
  
  // Probability of mild reported 
   vector[J] Pi; 
   // real sd_p = exp(-0.5 * psi_p);
   vector[J] logit_pi;
   
  logit_pi[1] = sd_p * z_logit_pi[1];
  for (j in 2:J){
    logit_pi[j] =  logit_pi[j - 1] + sd_p * z_logit_pi[j];}
   
  // Probability of severe reported
   matrix[I, J] tau;
   // real<lower=0> sd_tau = exp(-0.5*psi_tau);
   matrix[I, J] logit_tau; 
  
  for (i in 1:I){
      logit_tau[i,1] = tau_mean[i] + sd_tau * z_logit_tau[i,1];
  for (j in 2:J){
    logit_tau[i,j] = tau_mean[i] + (logit_tau[i,j - 1]-tau_mean[i])*rho_tau + sd_tau * z_logit_tau[i,j];}}

  array[I, J] simplex[3] p;
for (j in 1:J){ 
    Pi[j] = inv_logit(logit_pi[j]);
  for (i in 1:I){
    tau[i,j] = inv_logit(logit_tau[i, j]);
    p[i, j][1] = tau[i,j];
    p[i, j][2] =  (1-tau[i,j])*Pi[j];
    p[i, j][3] = (1-tau[i,j])*(1-Pi[j]);
  }}
    
}

model {
    // Prior on upper triangle
  // to_vector(C_raw) ~ normal(1.5,1);  
  // Prior on logit(pi)
  z_logit_pi ~ normal(0, 1);
  to_vector(z_logit_tau) ~ normal(0,1);
  
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
      target += multinomial_log_approx(A[i, j], Y[i, j], X[i, j], p[i,j]);
    }
  }
  
  
  // priors
  for (i in 1:I)
    C[i] ~ dirichlet(rep_vector(1.0, I));  // flat prior 
      
  // target += log(0.5 * phi_p) - phi_p*exp(-0.5*psi_p) - 0.5*psi_p;
  sd_p ~ exponential(phi_p); 
  sd_tau ~ exponential(phi_tau); 
  //admissions 
  // target += log(0.5 * phi_tau) - phi_tau*exp(-0.5*psi_tau) - 0.5*psi_tau;
  // rho_pi ~ beta(2,2);         //prior for geometric delay parameter
  rho_tau ~ beta(1,3);         //prior for geometric delay parameter
  tau_mean ~ normal(tau_mean_prior, 0.5);
}
