################################################################################
# Code to test various functions of beastt - normal endpoint
################################################################################

### Stan model - sigma2 known, mixture prior with a t component (e.g., power prior) and a normal component
BDB_stan_sigma2_known <- "data {
  int<lower = 0> N;                     // sample size of the internal control arm
  vector[N] y;                          // N x 1 vector of normal responses
  int<lower = 0> nu_pp;                 // degrees of freedom for t power prior
  real theta_pp;                        // mean hyperparameter for t power prior
  real<lower = 0> tau_pp;               // scale hyperparameter for t power prior
  real theta_v;                         // mean hyperparameter for normal vague prior
  real<lower = 0> tau_v;                // standard deviation hyperparameter for normal vague prior
  real<lower = 0, upper = 1> w;         // prior weight associated with power prior of RMP
  real<lower = 0> sigma_IC;             // standard deviation of internal control arm
}
parameters {
  real muC;                             // mean of control arm
}
model {
  // Prior distributions
  target += log_mix(w, student_t_lpdf(muC|nu_pp, theta_pp, tau_pp), normal_lpdf(muC|theta_v, tau_v));   // RMP with K=2 componenents
  target += -2 * log(sigma_IC);        // improper prior for internal control standard deviation

  // Likelihood
  y ~ normal(muC, sigma_IC);           // specify standard deviation in normal function
}
"


### Stan model - sigma2 unknown, mixture prior with a t component (e.g., power prior) and a normal component
BDB_stan_sigma2_unknown <- "data {
  int<lower = 0> N;                     // sample size of the internal control arm
  vector[N] y;                          // N x 1 vector of normal responses
  int<lower = 0> nu_pp;                 // degrees of freedom for t power prior
  real theta_pp;                        // mean hyperparameter for t power prior
  real<lower = 0> tau_pp;               // scale hyperparameter for t power prior
  real theta_v;                         // mean hyperparameter for normal vague prior
  real<lower = 0> tau_v;                // standard deviation hyperparameter for normal vague prior
  real<lower = 0, upper = 1> w;         // prior weight associated with power prior of RMP
}
parameters {
  real muC;                             // mean of control arm
  real<lower = 0> sigma_IC;             // standard deviation of internal control arm
}
model {
  // Prior distributions
  target += log_mix(w, student_t_lpdf(muC|nu_pp, theta_pp, tau_pp), normal_lpdf(muC|theta_v, tau_v));   // RMP with K=2 componenents
  target += -2 * log(sigma_IC);        // improper prior for internal control standard deviation

  // Likelihood
  y ~ normal(muC, sigma_IC);           // specify standard deviation in normal function
}
"
