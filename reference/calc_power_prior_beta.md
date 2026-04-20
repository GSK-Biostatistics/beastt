# Calculate Power Prior Beta

Calculate a (potentially inverse probability weighted) beta power prior
for the control response rate using external control data.

## Usage

``` r
calc_power_prior_beta(external_data, response, prior)
```

## Arguments

- external_data:

  This can either be a `prop_scr_obj` created by calling
  `create_prop_scr()` or a tibble of the external data. If it is just a
  tibble the weights will be assumed to be 1.

- response:

  Name of response variable

- prior:

  A beta distributional object that is the initial prior for the control
  response rate before the external control data are observed

## Value

Beta power prior object

## Details

Weighted participant-level response data from the control arm of an
external study are incorporated into an inverse probability weighted
(IPW) power prior for the control response rate \\\theta_C\\. When
borrowing information from an external control arm of size \\N\_{EC}\\,
the components of the IPW power prior for \\\theta_C\\ are defined as
follows:

- Initial prior::

  \$\$\theta_C \sim \mbox{Beta}(\nu_0, \phi_0)\$\$

- IPW likelihood of the external response data \\\boldsymbol{y}\_E\\
  with weights \\\hat{\boldsymbol{a}}\_0\\::

  \$\$\mathcal{L}\_E(\theta_C \mid \boldsymbol{y}\_E,
  \hat{\boldsymbol{a}}\_0) \propto \exp \left( \sum\_{i=1}^{N\_{EC}}
  \hat{a}\_{0i} \left\[ y_i \log(\theta_C) + (1 - y_i) \log(1 -
  \theta_C) \right\] \right)\$\$

- IPW power prior::

  \$\$\theta_C \mid \boldsymbol{y}\_E, \hat{\boldsymbol{a}}\_0 \sim
  \mbox{Beta} \left( \sum\_{i=1}^{N\_{EC}} \hat{a}\_{0i} y_i + \nu_0,
  \sum\_{i=1}^{N\_{EC}} \hat{a}\_{0i} (1 - y_i) + \phi_0 \right)\$\$

Defining the weights \\\hat{\boldsymbol{a}}\_0\\ to equal 1 results in a
conventional beta power prior.

## See also

Other power prior:
[`calc_power_prior_norm()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_norm.md),
[`calc_power_prior_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_weibull.md)

## Examples

``` r
library(distributional)
library(dplyr)
# This function can be used directly on the data
calc_power_prior_beta(external_data = ex_binary_df,
                      response = y,
                      prior = dist_beta(0.5, 0.5))
#> <distribution[1]>
#> [1] Beta(78, 74)

# Or this function can be used with a propensity score object
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                        external_df = ex_binary_df,
                        id_col = subjid,
                        model = ~ cov1 + cov2 + cov3 + cov4)

calc_power_prior_beta(ps_obj,
                      response = y,
                      prior = dist_beta(0.5, 0.5))
#> <distribution[1]>
#> [1] Beta(45, 37)
```
