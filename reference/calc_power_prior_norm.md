# Calculate Power Prior Normal

Calculate a (potentially inverse probability weighted) normal power
prior using external data.

## Usage

``` r
calc_power_prior_norm(
  external_data,
  response,
  prior = NULL,
  external_sd = NULL
)
```

## Arguments

- external_data:

  This can either be a `prop_scr_obj` created by calling
  `create_prop_scr()` or a tibble of the external data. If it is just a
  tibble the weights will be assumed to be 1. Only the external data for
  the arm(s) of interest should be included in this object (e.g.,
  external control data if creating a power prior for the control mean)

- response:

  Name of response variable

- prior:

  Either `NULL` or a normal distributional object that is the initial
  prior for the parameter of interest (e.g., control mean) before the
  external data are observed

- external_sd:

  Standard deviation of external response data if assumed known. It can
  be left as `NULL` if assumed unknown

## Value

Normal power prior object

## Details

Weighted participant-level response data from an external study are
incorporated into an inverse probability weighted (IPW) power prior for
the parameter of interest \\\theta\\ (e.g., the control mean if
borrowing from an external control arm). When borrowing information from
an external dataset of size \\N\_{E}\\, the IPW likelihood of the
external response data \\\boldsymbol{y}\_E\\ with weights
\\\hat{\boldsymbol{a}}\_0\\ is defined as

\$\$\mathcal{L}\_E(\theta \mid \boldsymbol{y}\_E,
\hat{\boldsymbol{a}}\_0, \sigma\_{E}^2) \propto \exp \left( -\frac{1}{2
\sigma\_{E}^2} \sum\_{i=1}^{N\_{E}} \hat{a}\_{0i} (y_i - \theta)^2
\right).\$\$

The `prior` argument should be either a distributional object with a
family type of `normal` or `NULL`, corresponding to the use of a normal
initial prior or an improper uniform initial prior (i.e., \\\pi(\theta)
\propto 1\\), respectively.

The `external_sd` argument can be a positive value if the external
standard deviation is assumed known or left as `NULL` otherwise. If
`external_sd = NULL`, then `prior` must be `NULL` to indicate the use of
an improper uniform initial prior for \\\theta\\, and an improper prior
is defined for the unknown external standard deviation such that
\\\pi(\sigma_E^2) \propto (\sigma_E^2)^{-1}\\. The details of the IPW
power prior for each case are as follows:

- `external_sd = positive value` (\\\sigma_E^2\\ known)::

  With either a proper normal or an improper uniform initial prior, the
  IPW weighted power prior for \\\theta\\ is a normal distribution.

- `external_sd = NULL` (\\\sigma_E^2\\ unknown)::

  With improper priors for both \\\theta\\ and \\\sigma_E^2\\, the
  marginal IPW weighted power prior for \\\theta\\ after integrating
  over \\\sigma_E^2\\ is a non-standardized \\t\\ distribution.

Defining the weights \\\hat{\boldsymbol{a}}\_0\\ to equal 1 results in a
conventional normal (or \\t\\) power prior if the external standard
deviation is known (unknown).

## See also

Other power prior:
[`calc_power_prior_beta()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_beta.md),
[`calc_power_prior_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_weibull.md)

## Examples

``` r
library(distributional)
library(dplyr)
# This function can be used directly on the data
calc_power_prior_norm(ex_norm_df,
                      response = y,
                      prior = dist_normal(0.5, 10),
                      external_sd = 0.15)
#> <distribution[1]>
#> [1] N(0.59, 0.00015)

# Or this function can be used with a propensity score object
ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                        external_df = ex_norm_df,
                        id_col = subjid,
                        model = ~ cov1 + cov2 + cov3 + cov4)
calc_power_prior_norm(ps_obj,
                      response = y,
                      prior = dist_normal(0.5, 10),
                      external_sd = 0.15)
#> <distribution[1]>
#> [1] N(0.69, 0.00037)
```
