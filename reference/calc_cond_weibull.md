# Calculate Conditional Drift and Treatment Effect for Time-to-Event Outcome Models

In order to properly generate time-to-event (TTE) outcome data for the
internal trial as part of a simulation study that investigates inverse
probability weighting, we need to translate the desired marginal drift
and treatment effect to the corresponding conditional drift and
treatment effect that can then be added into a TTE outcome model (e.g.,
Weibull proportional hazards regression model) used to simulate response
data.

## Usage

``` r
calc_cond_weibull(
  population,
  weibull_ph_mod,
  marg_drift,
  marg_trt_eff,
  analysis_time
)
```

## Arguments

- population:

  A very large data frame (e.g., number of rows \\\ge\\ 100,000) where
  the columns correspond to the covariates defined in the `survreg`
  object for the Weibull proportional hazards model. This data frame
  should be constructed to represent the population of the internal
  trial according to the assumed covariate distributions (possibly
  imbalanced from the external data).

- weibull_ph_mod:

  `survreg` object corresponding to a Weibull proportional hazards model
  fit using the external data

- marg_drift:

  Vector of marginal drift values

- marg_trt_eff:

  Vector of marginal treatment effect values

- analysis_time:

  A single time point when survival probabilities will be calculated

## Value

tibble of all combinations of the marginal drift and treatment effect.
For each row the conditional drift and treatment effect has been
calculated as well as the true marginal survival probabilities at time
`t` for the control and treatment populations.

## Details

In simulation studies that investigate the properties of inverse
probability weighted Bayesian dynamic borrowing, scenarios should be
considered in which the underlying survival probabilities at some
prespecified time \\t\\ (`analysis_time`) for the internal and external
control populations differ by varying amounts due to unmeasured
confounding (i.e., drift, where positive values indicate a higher
survival probability for the internal population). While values of drift
and treatment effect (i.e., difference between the survival
probabilities at time \\t\\ for the treated and control populations) can
be defined on the marginal scale for simulation studies, we must first
convert these values to the conditional scale and then include these
terms, along with covariates, in a Weibull proportional hazards (PH)
regression outcome model when generating time-to-event (TTE) data for
the internal arms. Doing so allows us to assume a relationship between
the covariates and the response variable while properly accounting for
drift and treatment effect.

To identify the conditional drift and treatment effect that correspond
to specified values of marginal drift and treatment effect, we first
bootstrap covariate vectors from the external data (e.g., \\N \ge
100,000\\) to construct a "population" that represents both the internal
trial (possibly incorporating intentional covariate imbalance) and the
external trial *after* standardizing it to match the covariate
distributions of the internal trial (allowing us to control for measured
confounding from potential imbalance in the covariate distributions).
Measured confounding can be incorporated into the data generation by
bootstrapping a very large data frame (`population`) in which the
distribution of at least one covariate is intentionally varied from that
of the external data; additional *unmeasured* drift can be incorporated
through the translation of specified marginal values (`marg_drift`) to
conditional values.

Let \\\Delta\\ and \\\delta\\ denote the marginal and conditional drift,
respectively. For a specified value of \\\Delta\\, we can identify the
corresponding \\\delta\\ as the value that, when added as an additional
term in the Weibull PH model survival function (i.e., additive change in
the intercept) for each individual in the population,
increases/decreases the population-averaged conditional probabilities of
survival at time \\t\\ by an amount approximately equal to \\\Delta\\.
That is, the optimal \\\delta\\ minimizes

\$\$\left\| \left( \frac{1}{N} \sum\_{i=1}^N \exp \left( -\left\\ \exp
\left( \boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta
\right) \times t \right\\^{\alpha\_{EC}} \right) - \frac{1}{N}
\sum\_{i=1}^N \exp \left( -\left\\ \exp \left( \boldsymbol{x}\_i^\prime
\boldsymbol{\beta}\_{EC} \right) \times t \right\\^{\alpha\_{EC}}
\right) \right) - \Delta \right\|,\$\$

where \\\alpha\_{EC}\\ is the Weibull shape parameter,
\\\boldsymbol{\beta}\_{EC}\\ is a vector of regression coefficients, and
\\\boldsymbol{x}\_i\\ is a vector of covariates (including an intercept
term) from the bootstrapped population of size \\N\\. We note that
\\\alpha\_{EC} = 1/\sigma\_{EC}\\ and \\\boldsymbol{\beta}\_{EC} =
-\boldsymbol{\xi}\_{EC}\\ are calculated as functions of the scale
parameter (\\\sigma\_{EC}\\) and coefficients
(\\\boldsymbol{\xi}\_{EC}\\) estimated by the `survreg` object that was
fit to the external data, and we assume here that these estimates are
the "true" shape and covariate effects when generating response data. In
the formula above, the first and second terms correspond to the
population-averaged conditional survival functions (i.e., the marginal
survival probabilities) at time \\t\\ for the internal control
population with drift and the external control population (with
covariate distributions standardized to match the internal trial),
respectively.

If we now denote the marginal and conditional treatment effect by
\\\Gamma\\ and \\\gamma\\, respectively, we can use a similar process to
identify the optimal \\\gamma\\ that approximately corresponds to the
specified value of \\\Gamma\\, which is done by minimizing the
following:

\$\$\left\| \left( \frac{1}{N} \sum\_{i=1}^N \exp \left( -\left\\ \exp
\left( \boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta +
\gamma \right) \times t \right\\^{\alpha\_{EC}} \right) - \frac{1}{N}
\sum\_{i=1}^N \exp \left( -\left\\ \exp \left( \boldsymbol{x}\_i^\prime
\boldsymbol{\beta}\_{EC} + \delta \right) \times t
\right\\^{\alpha\_{EC}} \right) \right) - \Gamma \right\|,\$\$

where the first term is the average of the conditional survival
functions (i.e., the marginal survival probabilities) at time \\t\\ for
the internal treated population.

See
[here](https://github.com/GSK-Biostatistics/beastt/blob/e2b41fe90f639924d10c0d94ceff04a74d0ce617/inst/templates/tte-template.R)
for a simulation example with a time-to-event outcome.

## Examples

``` r
library(dplyr)
library(survival)
# Model "true" regression coefficients using the external data
weibull_ph_mod <- survreg(Surv(y, event) ~ cov1 + cov2 + cov3 + cov4, data = ex_tte_df,
                          dist = "weibull")

# Bootstrap internal control "population" with imbalance w.r.t. covariate 2
pop_int_ctrl <- bootstrap_cov(ex_tte_df, n = 100000, imbal_var = cov2,
                              imbal_prop = 0.25, ref_val = 0) |>
  select(c(cov1, cov2, cov3, cov4))     # keep only covariate columns

# Convert the marginal drift and treatment effects to conditional
calc_cond_weibull(population = pop_int_ctrl, weibull_ph_mod,
                  marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .10),
                  analysis_time = 12)
#> # A tibble: 6 × 6
#>   marg_drift marg_trt_eff conditional_drift true_control_surv_prob
#>        <dbl>        <dbl>             <dbl>                  <dbl>
#> 1       -0.1          0               0.210                  0.367
#> 2       -0.1          0.1             0.210                  0.367
#> 3        0            0               0                      0.467
#> 4        0            0.1             0                      0.467
#> 5        0.1          0              -0.223                  0.567
#> 6        0.1          0.1            -0.223                  0.567
#> # ℹ 2 more variables: conditional_trt_eff <dbl>, true_trt_surv_prob <dbl>
```
