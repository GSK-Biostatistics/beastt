data {
  int<lower=0> N;         // sample size
  vector[N] y;            // observed time (event or censored)
  vector[N] e;            // event indicator (1: event; 0: censored)

  int<lower=0> K;         // number of times to compute survival probabilities
  vector[K]    times;     // vector of survival probability times

  // Robust mixture prior hyperparameters
  matrix[2,2] cov_inf;        // covariance matrix of informative prior component
  matrix[2,2] cov_vague;      // covariance matrix of vague prior component

  vector[2] mean_inf;           // mean vector of informative prior component
  vector[2] mean_vague;         // mean vector of vague prior component
  real<lower=0, upper=1> w0;    // prior mixture weight associated with informative component
}
parameters {
  real beta0;         // regression intercept parameter
  real log_alpha;     // log of the shape parameter
}
transformed parameters {
  real alpha = exp(log_alpha);    // shape parameter
}
model {
  vector[2] theta;
  vector[2] log_prior;
  theta[1] = log_alpha;
  theta[2] = beta0;

  // Robust mixture prior
  log_prior[1] = log(w0)   + multi_normal_lpdf(theta|mean_inf, cov_inf);
  log_prior[2] = log(1-w0) + multi_normal_lpdf(theta|mean_vague, cov_vague);
  target += log_sum_exp( log_prior );

  // Likelihood
  for (i in 1:N) {
    if (e[i] == 1) {
      target += weibull_lpdf(y[i] | alpha, exp(-beta0));
    } else {
      target += weibull_lccdf(y[i] | alpha, exp(-beta0));
    }
  }
}
generated quantities {
  vector[K] survProb;
  for (t in 1:K) {
    survProb[t] = exp(weibull_lccdf(times[t] | alpha, exp(-beta0)));
  }
}
