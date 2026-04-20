# Calculate Posterior Beta

Calculate a posterior distribution that is beta (or a mixture of beta
components). Only the relevant treatment arms from the internal dataset
should be read in (e.g., only the control arm if constructing a
posterior distribution for the control response rate).

## Usage

``` r
calc_post_beta(internal_data, response, prior)
```

## Arguments

- internal_data:

  A tibble of the internal data.

- response:

  Name of response variable

- prior:

  A distributional object corresponding to a beta distribution or a
  mixture distribution of beta components

## Value

distributional object

## Details

For a given arm of an internal trial (e.g., the control arm or an active
treatment arm) of size \\N_I\\, suppose the response data are binary
such that \\Y_i \sim \mbox{Bernoulli}(\theta)\\, \\i=1,\ldots,N_I\\. The
posterior distribution for \\\theta\\ is written as

\$\$\pi( \theta \mid \boldsymbol{y} ) \propto \mathcal{L}(\theta \mid
\boldsymbol{y}) \\ \pi(\theta),\$\$

where \\\mathcal{L}(\theta \mid \boldsymbol{y})\\ is the likelihood of
the response data from the internal arm and \\\pi(\theta)\\ is a prior
distribution on \\\theta\\ (either a beta distribution or a mixture
distribution with an arbitrary number of beta components). The posterior
distribution for \\\theta\\ is either a beta distribution or a mixture
of beta components depending on whether the prior is a single beta
distribution or a mixture distribution.

## Examples

``` r
library(dplyr)
library(distributional)
calc_post_beta(internal_data = filter(int_binary_df, trt == 1),
               response = y,
               prior = dist_beta(0.5, 0.5))
#> <distribution[1]>
#> [1] Beta(62, 20)
```
