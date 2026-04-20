# Simulate Event Times for Each Participant from a Weibull Proportional Hazards Regression Model

Simulate Event Times for Each Participant from a Weibull Proportional
Hazards Regression Model

## Usage

``` r
sim_weib_ph(weibull_ph_mod, samp_df, cond_drift = 0, cond_trt_effect = 0)
```

## Arguments

- weibull_ph_mod:

  `survreg` object corresponding to a Weibull proportional hazards model
  fit using the external data

- samp_df:

  Data frame of covariates corresponding to the sample arm (control or
  treated) for which event times should be simulated. The column names
  should correspond to the covariate names in the `survreg` object.

- cond_drift:

  Optional value of the conditional drift by which the intercept in the
  Weibull proportional hazards regression model should be
  increased/decreased to incorporate the impact of unmeasurable sources
  of drift. Default is 0.

- cond_trt_effect:

  Optional value of the conditional treatment effect by which the
  intercept in the Weibull proportional hazards regression model should
  be increased/decreased if simulating event data for a treated arm.
  Default is 0.

## Value

Vector of simulated event times from a Weibull proportional hazards
regression model

## Details

Simulate the event times for each participant using a Weibull
proportional hazards (PH) regression model. The "true" parameter values
for the Weibull shape \\\alpha\\ and the regression coefficients
\\\boldsymbol{\beta}\\ are assumed to be equal to the parameter
estimates from a `survreg` object (`weibull_ph_mod`) fit using external
data (note that the Weibull shape parameter \\\alpha\\ is defined as the
inverse of the scale parameter reported by `survreg`).

For participant \\i\\, let \\y_i\\ denote the time-to-event random
variable and \\\boldsymbol{x}\_i = \\x\_{i,1}, \ldots, x\_{i,p}\\\\ the
vector of \\p\\ covariates (row \\i\\ of `samp_df`) that correspond to
the \\(p+1)\\-dimensional vector of regression coefficients
\\\boldsymbol{\beta}\\. The density function of the Weibull PH
regression model is

\$\$f(y_i \mid \boldsymbol{x}\_i, \alpha, \boldsymbol{\beta}, \delta,
\gamma) = \left( \frac{\alpha}{\sigma_i} \right) \left(
\frac{y_i}{\sigma_i} \right)^{\alpha - 1} \exp \left( -\left(
\frac{y_i}{\sigma_i} \right)^\alpha \right),\$\$

where \\-\log(\sigma_i) = \beta_0 + \beta_1 x\_{i,1} + \ldots + \beta_p
x\_{i,p} + \delta + \gamma\\. Here, \\\delta\\ and \\\gamma\\ denote the
conditional drift (`cond_drift`) and conditional treatment effect
(`cond_trt_effect`), respectively, that can be calculated using
[`calc_cond_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_cond_weibull.md)
for desired values of the marginal drift and marginal treatment effect.

## Examples

``` r
library(dplyr)
library(survival)
# Model "true" regression coefficients and shape parameter using the external data
weibull_ph_mod <- survreg(Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = ex_tte_df,
                          dist = "weibull")

# Sample covariates for internal control arm via bootstrap from external data
samp_int_ctrl <- bootstrap_cov(ex_tte_df, n = 100) |>
  select(c(cov1, cov2, cov3, cov4))     # keep only covariate columns
tte_dat <- sim_weib_ph(weibull_ph_mod, samp_df = samp_int_ctrl)
```
