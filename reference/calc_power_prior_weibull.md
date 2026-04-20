# Calculate Power Prior Weibull

Calculate an approximate (potentially inverse probability weighted)
multivariate normal power prior for the log-shape and log-inverse-scale
parameters of a Weibull likelihood for external time-to-event control
data.

## Usage

``` r
calc_power_prior_weibull(
  external_data,
  response,
  event,
  intercept,
  shape,
  approximation = c("Laplace", "MCMC"),
  ...
)
```

## Arguments

- external_data:

  This can either be a `prop_scr_obj` created by calling
  `create_prop_scr()` or a tibble of the external data. If it is just a
  tibble the weights will be assumed to be 1. Only the external data for
  the arm(s) of interest should be included in this object (e.g.,
  external control data if creating a power prior for the control
  Weibull shape and intercept parameters)

- response:

  Name of response variable

- event:

  Name of event indicator variable (1: event; 0: censored)

- intercept:

  Normal distributional object that is the initial prior for the
  intercept (i.e., log-inverse-scale) parameter

- shape:

  Integer value that is the scale of the half-normal prior for the shape
  parameter

- approximation:

  Type of approximation to be used. Either `Laplace` or `MCMC`.
  `Laplace` is used by default because it is faster than `MCMC`.

- ...:

  Arguments passed to
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  (e.g. iter, chains).

## Value

Multivariate Normal Distributional Object

## Details

Weighted participant-level response data from the control arm of an
external study are incorporated into an approximated inverse probability
weighted (IPW) power prior for the parameter vector
\\\boldsymbol{\theta}\_C = \\\log(\alpha), \beta\\\\, where \\\beta =
-\log(\sigma)\\ is the intercept parameter of a Weibull proportional
hazards regression model and \\\alpha\\ and \\\sigma\\ are the Weibull
shape and scale parameters, respectively. When borrowing information
from an external dataset of size \\N\_{E}\\, the IPW likelihood of the
external response data \\\boldsymbol{y}\_E\\ with event indicators
\\\boldsymbol{\nu}\_E\\ and weights \\\hat{\boldsymbol{a}}\_0\\ is
defined as

\$\$\mathcal{L}\_E(\alpha, \sigma \mid \boldsymbol{y}\_E,
\boldsymbol{\nu}\_E, \hat{\boldsymbol{a}}\_0) \propto \prod\_{i=1}^{N_E}
\left\\ \left( \frac{\alpha}{\sigma} \right) \left( \frac{y_i}{\sigma}
\right)^{\alpha - 1} \exp \left( -\left( \frac{y_i}{\sigma}
\right)^\alpha \right) \right\\^{\hat{a}\_{0i} \nu_i} \left\\ \exp
\left( -\left( \frac{y_i}{\sigma} \right)^\alpha \right)
\right\\^{\hat{a}\_{0i}(1 - \nu_i)}.\$\$

The initial priors for the intercept parameter \\\beta\\ and the shape
parameter \\\alpha\\ are assumed to be normal and half-normal,
respectively, which can be defined using the `intercept` and `shape`
arguments.

The power prior for \\\boldsymbol{\theta}\_C\\ does not have a closed
form, and thus we approximate it via a bivariate normal distribution;
i.e., \$\$\boldsymbol{\theta}\_C \mid \boldsymbol{y}\_E,
\boldsymbol{\nu}\_E, \hat{\boldsymbol{a}}\_0 \\ \dot\sim \\ \mbox{MVN}
\left( \tilde{\boldsymbol{\mu}}\_0, \tilde{\boldsymbol{\Sigma}}\_0
\right).\$\$ If `approximation = Laplace`, then
\\\tilde{\boldsymbol{\mu}}\_0\\ is the mode vector of the IPW power
prior and \\\tilde{\boldsymbol{\Sigma}}\_0\\ is the negative inverse of
the Hessian of the log IPW power prior evaluated at the mode. If
`approximation = MCMC`, then MCMC samples are obtained from the IPW
power prior, and \\\tilde{\boldsymbol{\mu}}\_0\\ and
\\\tilde{\boldsymbol{\Sigma}}\_0\\ are the estimated mean vector and
covariance matrix of these MCMC samples. Note that the Laplace
approximation method is faster due to its use of optimization instead of
MCMC sampling.

The first element of the mean vector and the first row/column of
covariance matrix correspond to the log-shape parameter
(\\\log(\alpha)\\), and the second element corresponds to the intercept
(\\\beta\\, the log-inverse-scale) parameter.

## See also

Other power prior:
[`calc_power_prior_beta()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_beta.md),
[`calc_power_prior_norm()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_norm.md)

## Examples

``` r
library(distributional)
library(dplyr)
# This function can be used directly on the data
calc_power_prior_weibull(ex_tte_df,
                         response = y,
                         event = event,
                         intercept = dist_normal(0, 10),
                         shape = 50,
                         approximation = "Laplace")
#> <distribution[1]>
#> [1] MVN[2]

# Or this function can be used with a propensity score object
ps_obj <- calc_prop_scr(internal_df = filter(int_tte_df, trt == 0),
                        external_df = ex_tte_df,
                        id_col = subjid,
                        model = ~ cov1 + cov2 + cov3 + cov4)
calc_power_prior_weibull(ps_obj,
                         response = y,
                         event = event,
                         intercept = dist_normal(0, 10),
                         shape = 50,
                         approximation = "Laplace")
#> <distribution[1]>
#> [1] MVN[2]
```
