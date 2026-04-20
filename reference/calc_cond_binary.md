# Calculate Conditional Drift and Treatment Effect for Binary Outcome Models

In order to properly generate binary response data for the internal
trial as part of a simulation study that investigates inverse
probability weighting, we need to translate the desired marginal drift
and treatment effect to the corresponding conditional drift and
treatment effect that can then be added into a binary outcome model
(e.g., logistic regression model) used to simulate response data.

## Usage

``` r
calc_cond_binary(population, glm, marg_drift, marg_trt_eff)
```

## Arguments

- population:

  A very large data frame (e.g., number of rows \\\ge\\ 100,000) where
  the columns correspond to the covariates defined in the logistic
  regression model object. This data frame should be constructed to
  represent the population of the internal trial according to the
  assumed covariate distributions (possibly imbalanced from the external
  data).

- glm:

  Logistic regression model object fit using the external data

- marg_drift:

  Vector of marginal drift values

- marg_trt_eff:

  Vector of marginal treatment effect values

## Value

tibble of all combinations of the marginal drift and treatment effect.
For each row the conditional drift and treatment effect has been
calculated as well as the true response rates for the control and
treated populations.

## Details

In simulation studies that investigate the properties of inverse
probability weighted Bayesian dynamic borrowing, scenarios should be
considered in which the underlying response rates for the internal and
external control populations differ by varying amounts due to unmeasured
confounding (i.e., drift, where positive values indicate a higher
response rate for the internal population). While values of drift and
treatment effect (i.e., risk difference) can be defined on the marginal
scale for simulation studies, we must first convert these values to the
conditional scale and then include these terms, along with covariates,
in a logistic regression outcome model when generating response data for
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
term in the logistic regression model (i.e., change in the intercept)
for each individual in the population, increases/decreases the
population-averaged conditional probabilities of response by an amount
approximately equal to \\\Delta\\. That is, the optimal \\\delta\\
minimizes

\$\$\left\| \left( \frac{1}{N} \sum\_{i=1}^N \frac{\exp \left(
\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta \right)}{1 +
\exp\left(\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta
\right)} - \frac{1}{N} \sum\_{i=1}^N \frac{\exp \left(
\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} \right)}{1 + \exp
\left(\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} \right)}
\right) - \Delta \right\|,\$\$

where \\\boldsymbol{\beta}\_{EC}\\ is the vector of regression
coefficients from the logistic regression model (`glm`) fit to the
external control data (assumed here to be the "true" covariate effects
when generating response data) and \\\boldsymbol{x}\_i\\ is a vector of
covariates (including an intercept term) from the bootstrapped
population of size \\N\\. In the formula above, the first and second
terms correspond to the population-averaged conditional probabilities
(i.e., the marginal response rates) of the internal control population
with drift and the external control population (with covariate
distributions standardized to match the internal trial), respectively.

If we now denote the marginal and conditional treatment effect by
\\\Gamma\\ and \\\gamma\\, respectively, we can use a similar process to
identify the optimal \\\gamma\\ that approximately corresponds to the
specified value of \\\Gamma\\, which is done by minimizing the
following:

\$\$\left\| \left( \frac{1}{N} \sum\_{i=1}^N \frac{\exp \left(
\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta + \gamma
\right)}{1 + \exp\left(\boldsymbol{x}\_i^\prime
\boldsymbol{\beta}\_{EC} + \delta + \gamma \right)} - \frac{1}{N}
\sum\_{i=1}^N \frac{\exp \left( \boldsymbol{x}\_i^\prime
\boldsymbol{\beta}\_{EC} + \delta \right)}{1 + \exp
\left(\boldsymbol{x}\_i^\prime \boldsymbol{\beta}\_{EC} + \delta
\right)} \right) - \Gamma \right\|,\$\$

where the first term is the average of the conditional probabilities of
response (i.e., the marginal response rate) of the internal treated
population.

See
[here](https://github.com/GSK-Biostatistics/beastt/blob/e2b41fe90f639924d10c0d94ceff04a74d0ce617/inst/templates/binary_template.R)
for a simulation example with a binary outcome.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
# Model "true" regression coefficients using the external data
logit_mod <- glm(y ~ cov1 + cov2 + cov3 + cov4, data = ex_binary_df, family = binomial)

# Bootstrap internal control "population" with imbalance w.r.t. covariate 2
pop_int_ctrl <- bootstrap_cov(ex_binary_df, n = 100000, imbal_var = cov2,
                              imbal_prop = 0.25, ref_val = 0) |>
                              select(-subjid, -y)  # keep only covariate columns

# Convert the marginal drift and treatment effects to conditional
calc_cond_binary(population = pop_int_ctrl, glm = logit_mod,
                 marg_drift = c(-.1, 0, .1), marg_trt_eff = c(0, .15))
#> # A tibble: 6 × 6
#>   marg_drift marg_trt_eff conditional_drift true_control_RR conditional_trt_eff
#>        <dbl>        <dbl>             <dbl>           <dbl>               <dbl>
#> 1       -0.1         0               -0.413           0.446               0    
#> 2       -0.1         0.15            -0.413           0.446               0.622
#> 3        0           0                0               0.546               0    
#> 4        0           0.15             0               0.546               0.660
#> 5        0.1         0                0.428           0.646               0    
#> 6        0.1         0.15             0.428           0.646               0.776
#> # ℹ 1 more variable: true_trt_RR <dbl>
```
