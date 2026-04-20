# Robustify Multivariate Normal Distributions

Adds vague normal component, where the level of vagueness is controlled
by the `n` parameter

## Usage

``` r
robustify_mvnorm(prior, n, weights = c(0.5, 0.5))
```

## Arguments

- prior:

  Multivariate Normal distributional object

- n:

  Number of theoretical participants (or events, for time-to-event data)

- weights:

  Vector of weights, where the first number corresponds to the
  informative component and the second is the vague

## Value

mixture distribution

## Details

In cases with a time-to-event endpoint, a robust mixture prior can be
created by adding a vague multivariate normal component to any
multivariate normal prior with mean vector \\\boldsymbol{\mu}\\ and
covariance matrix \\\boldsymbol{\Sigma}\\. The vague component is
calculated to have the same mean vector \\\boldsymbol{\mu}\\ and
covariance matrix equal to \\\boldsymbol{\Sigma} \times n\\, where `n`
is the specified number of theoretical events.

## Examples

``` r
library(distributional)
robustify_mvnorm(
      dist_multivariate_normal(mu = list(c(1, 0)), sigma = list(c(10, 5))),
       n = 15)
#> <distribution[1]>
#> [1] mixture(0.5*MVN[2], 0.5*MVN[2])
```
