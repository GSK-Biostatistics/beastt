# Calculate Posterior Normal

Calculate a posterior distribution that is normal (or a mixture of
normal components). Only the relevant treatment arms from the internal
dataset should be read in (e.g., only the control arm if constructing a
posterior distribution for the control mean).

## Usage

``` r
calc_post_norm(internal_data, response, prior, internal_sd = NULL)
```

## Arguments

- internal_data:

  A tibble of the internal data.

- response:

  Name of response variable

- prior:

  A distributional object corresponding to a normal distribution, a t
  distribution, or a mixture distribution of normal and/or t components

- internal_sd:

  Standard deviation of internal response data if assumed known. It can
  be left as `NULL` if assumed unknown

## Value

distributional object

## Details

For a given arm of an internal trial (e.g., the control arm or an active
treatment arm) of size \\N_I\\, suppose the response data are normally
distributed such that \\Y_i \sim N(\theta, \sigma_I^2)\\,
\\i=1,\ldots,N_I\\. If \\\sigma_I^2\\ is assumed known, the posterior
distribution for \\\theta\\ is written as

\$\$\pi( \theta \mid \boldsymbol{y}, \sigma\_{I}^2 ) \propto
\mathcal{L}(\theta \mid \boldsymbol{y}, \sigma\_{I}^2) \\
\pi(\theta),\$\$

where \\\mathcal{L}(\theta \mid \boldsymbol{y}, \sigma\_{I}^2)\\ is the
likelihood of the response data from the internal arm and
\\\pi(\theta)\\ is a prior distribution on \\\theta\\ (either a normal
distribution, a \\t\\ distribution, or a mixture distribution with an
arbitrary number of normal and/or \\t\\ components). Any \\t\\
components of the prior for \\\theta\\ are approximated with a mixture
of two normal distributions.

If \\\sigma_I^2\\ is unknown, the marginal posterior distribution for
\\\theta\\ is instead written as

\$\$\pi( \theta \mid \boldsymbol{y} ) \propto \left\\ \int_0^\infty
\mathcal{L}(\theta, \sigma\_{I}^2 \mid \boldsymbol{y}) \\
\pi(\sigma\_{I}^2) \\ d\sigma\_{I}^2 \right\\ \times \pi(\theta).\$\$

In this case, the prior for \\\sigma_I^2\\ is chosen to be
\\\pi(\sigma\_{I}^2) = (\sigma_I^2)^{-1}\\ such that \\\left\\
\int_0^\infty \mathcal{L}(\theta, \sigma\_{I}^2 \mid \boldsymbol{y}) \\
\pi(\sigma\_{I}^2) \\ d\sigma\_{I}^2 \right\\\\ becomes a
non-standardized \\t\\ distribution. This integrated likelihood is then
approximated with a mixture of two normal distributions.

If `internal_sd` is supplied a positive value and `prior` corresponds to
a single normal distribution, then the posterior distribution for
\\\theta\\ is a normal distribution. If `internal_sd = NULL` or if other
types of prior distributions are specified (e.g., mixture or t
distribution), then the posterior distribution is a mixture of normal
distributions.

## Examples

``` r
library(distributional)
library(dplyr)
post_treated <- calc_post_norm(internal_data = filter(int_norm_df, trt == 1),
                               response = y,
                               prior = dist_normal(0.5, 10),
                               internal_sd = 0.15)
```
