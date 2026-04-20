# Calculate Posterior Weibull

Calculate a posterior distribution for the probability of surviving past
a given analysis time(s) for time-to-event data with a Weibull
likelihood. Only the relevant treatment arms from the internal dataset
should be read in (e.g., only the control arm if constructing a
posterior distribution for the control survival probability).

## Usage

``` r
calc_post_weibull(internal_data, response, event, prior, analysis_time, ...)
```

## Arguments

- internal_data:

  This can either be a propensity score object or a tibble of the
  internal data.

- response:

  Name of response variable

- event:

  Name of event indicator variable (1: event; 0: censored)

- prior:

  A distributional object corresponding to a multivariate normal
  distribution or a mixture of 2 multivariate normals. The first element
  of the mean vector and the first row/column of covariance matrix
  correspond to the log-shape parameter, and the second element
  corresponds to the intercept (i.e., log-inverse-scale) parameter.

- analysis_time:

  Vector of time(s) when survival probabilities will be calculated

- ...:

  rstan sampling option. This overrides any default options. For more
  information, see
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

## Value

stan posterior object

## Details

For a given arm of an internal trial (e.g., the control arm or an active
treatment arm) of size \\N_I\\, suppose the response data are time to
event such that \\Y_i \sim \mbox{Weibull}(\alpha, \sigma)\\, where

\$\$f(y_i \mid \alpha, \sigma) = \left( \frac{\alpha}{\sigma} \right)
\left( \frac{y_i}{\sigma} \right)^{\alpha - 1} \exp \left( -\left(
\frac{y_i}{\sigma} \right)^\alpha \right),\$\$

\\i=1,\ldots,N_I\\. Define \\\boldsymbol{\theta} = \\\log(\alpha),
\beta\\\\ where \\\beta = -\log(\sigma)\\ is the intercept parameter of
a Weibull proportional hazards regression model. The posterior
distribution for \\\boldsymbol{\theta}\\ is written as

\$\$\pi( \boldsymbol{\theta} \mid \boldsymbol{y}, \boldsymbol{\nu} )
\propto \mathcal{L}(\boldsymbol{\theta} \mid \boldsymbol{y},
\boldsymbol{\nu}) \\ \pi(\boldsymbol{\theta}),\$\$

where \\\mathcal{L}(\boldsymbol{\theta} \mid \boldsymbol{y},
\boldsymbol{\nu}) = \prod\_{i=1}^{N_I} f(y_i \mid
\boldsymbol{\theta})^{\nu_i} S(y_i \mid \boldsymbol{\theta})^{1 -
\nu_i}\\ is the likelihood of the response data from the internal arm
with event indicator \\\boldsymbol{\nu}\\ and survival function \\S(y_i
\mid \boldsymbol{\theta}) = 1 - F(y_i \mid \boldsymbol{\theta})\\, and
\\\pi(\boldsymbol{\theta})\\ is a prior distribution on
\\\boldsymbol{\theta}\\ (either a multivariate normal distribution or a
mixture of two multivariate normal distributions). Note that the
posterior distribution for \\\boldsymbol{\theta}\\ does not have a
closed form, and thus MCMC samples for \\\log(\alpha)\\ and \\\beta\\
are drawn from the posterior. These MCMC samples are used to construct
samples from the posterior distribution for the probability of surviving
past the specified analysis time(s) for the given arm.

To construct a posterior distribution for the treatment difference
(i.e., the difference in survival probabilities at the specified
analysis time), obtain MCMC samples from the posterior distributions for
the survival probabilities under each arm and then subtract the
control-arm samples from the treatment-arm samples.

## Examples

``` r
library(distributional)
library(dplyr)
library(rstan)
#> Loading required package: StanHeaders
#> 
#> rstan version 2.32.7 (Stan version 2.32.2)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
#> For within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,
#> change `threads_per_chain` option:
#> rstan_options(threads_per_chain = 1)
mvn_prior <- dist_multivariate_normal(
   mu = list(c(0.3, -2.6)),
   sigma = list(matrix(c(1.5, 0.3, 0.3, 1.1), nrow = 2)))
post_treated <- calc_post_weibull(filter(int_tte_df, trt == 1),
                                  response = y,
                                  event = event,
                                  prior = mvn_prior,
                                  analysis_time = 12,
                                  warmup = 5000,
                                  iter = 15000)
# Extract MCMC samples of survival probabilities at time t=12
surv_prob_treated <- as.data.frame(extract(post_treated,
                                   pars = c("survProb")))[,1]
```
