################################################################################
# Code to test various functions of beastt - TTE endpoint
################################################################################

### Stan model - power prior
BDB_stan_pp <- "data {
  int<lower=0> N;           // sample size
  vector[N] y;              // observation time
  vector[N] e;              // event indicator
  vector[N] wgt;            // weight

  real          beta_mean;  // mean for normal prior on the regression parameters;
  real<lower=0> beta_sd;    // standard deviation for normal prior on the regression parameters;

  real          shape_mean; // mean of half normal prior on weibull shape parameter
  real<lower=0> shape_sd;   // sd of half normal prior on weibull shape parameter
}

parameters {
  real beta;                // intercept parameter
  real<lower=0> shape;      // shape parameter
}
model {

  //------------------------
  // Priors
  //------------------------
  shape     ~ normal(shape_mean, shape_sd);
  beta      ~ normal(beta_mean,  beta_sd);

  for (i in 1:N)
  {
    if (e[i] == 1) {
      target += wgt[i]*weibull_lpdf(y[i] | shape, exp(-beta));
    } else {
      target += wgt[i]*weibull_lccdf(y[i] | shape, exp(-beta));
    }
  }
}
generated quantities
{

}
"


### Stan model - posterior
BDB_stan_post <- "data {
  int<lower=0> N;         // sample size
  vector[N] y;            // observation time
  vector[N] e;            // event indicator

  int<lower=0> T;         // Number of times to compute survival probabilities
  vector[T]    times;     // vector of survival probability times

  // robust mixture prior hyperparameters

  matrix[2,2] cov_inf;        // informative prior component covariance matrix
  matrix[2,2] cov_vague;      // vague prior component covariance matrix

  vector[2] mu_inf;           // informative prior component mean vector
  vector[2] mu_vague;         // vague prior component mean vector
  real<lower=0,upper=1> w0;   // prior mixture weight

}

parameters {
  real beta;                                // intercept parameter
  real logShape;                                // shape parameter
}

transformed parameters {
 real shape = exp(logShape);
}

model {

  vector[2] theta;
  vector[2] logPrior;
  //------------------------
  // Robust Mixture Prior
  //------------------------
  theta[1] = beta;
  theta[2] = logShape;


  logPrior[1] = log(w0)   + multi_normal_lpdf(theta|mu_inf,  cov_inf);
  logPrior[2] = log(1-w0) + multi_normal_lpdf(theta|mu_vague,cov_vague);
  target += log_sum_exp( logPrior );

  for (i in 1:N)
  {
    if (e[i] == 1) {
      target += weibull_lpdf(y[i] | shape, exp(-beta));
    } else {
      target += weibull_lccdf(y[i] | shape, exp(-beta));
    }
  }
}
generated quantities
{
  vector[T] survProb;
  for (t in 1:T)
  {
    survProb[t] = exp(weibull_lccdf( times[t] | shape, exp(-beta)));
  }
}
"
