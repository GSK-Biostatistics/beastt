# Robustify Normal Distributions

Adds vague normal component, where the level of vagueness is controlled
by the `n` parameter

## Usage

``` r
robustify_norm(prior, n, weights = c(0.5, 0.5))
```

## Arguments

- prior:

  Normal or Multivariate Normal distributional object

- n:

  Number of theoretical participants (or events, for time-to-event data)

- weights:

  Vector of weights, where the first number corresponds to the
  informative component and the second is the vague

## Value

mixture distribution

## Details

In cases with a normal endpoint, a robust mixture prior can be created
by adding a vague normal component to any normal prior with mean
\\\theta\\ and variance \\\sigma^2\\.The vague component is calculated
to have the same mean \\\theta\\ and variance equal to \\\sigma^2 \times
n\\, where `n` is the specified number of theoretical participants. If
robustifying a normal power prior that was calculated from external
control data and `n` is defined as the number of external control
participants, and the vague component would then correspond to one
external control participant's worth of data.

## Examples

``` r
library(distributional)
robustify_norm(dist_normal(0,1), n = 15)
#> <distribution[1]>
#> [1] mixture(0.5*N(0, 1), 0.5*N(0, 15))
```
