//
// Weibull Power Prior
//

data {
  int<lower=0> N;            // sample size
  vector[N] y;               // observed time (event or censored)
  vector[N] e;               // event indicator (1: event; 0: censored)
  vector[N] wgt;             // weight (e.g., inverse probability weight)

  real          beta0_mean;  // mean of normal prior on the regression intercept parameter
  real<lower=0> beta0_sd;    // standard deviation of normal prior on the regression intercept parameter

  real<lower=0> alpha_scale; // scale of half-normal prior on Weibull shape parameter
}
parameters {
  real beta0;               // regression intercept parameter
  real<lower=0> alpha;      // shape parameter
}
model {
  // Prior distributions
  target += normal_lpdf(alpha | 0, alpha_scale);
  target += normal_lpdf(beta0 | beta0_mean, beta0_sd);

  // Weighted likelihood
  for (i in 1:N) {
    if (e[i] == 1) {
      target += wgt[i]*weibull_lpdf(y[i] | alpha, exp(-beta0));
    } else {
      target += wgt[i]*weibull_lccdf(y[i] | alpha, exp(-beta0));
    }
  }
}
