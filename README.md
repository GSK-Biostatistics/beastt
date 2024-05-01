
# beastt

## Bayesian Evaluation, Analysis, and Simulation Software Tools for Trials (BEASTT)

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/beastt)](https://CRAN.R-project.org/package=beastt)
[![R-CMD-check](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

Welcome to the “beastt” package! This R package is designed to assist
users in performing Bayesian Dynamic Borrowing with covariate adjustment
for simulations and data analyses in clinical trials.

## Installation

You can install the development version of beastt from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GSK-Biostatistics/beastt")
```

## Functions

### `calc_prop_scr()`

Calculate propensity scores using this function and obtain the inverse
probability weight calculated for each subject. Propensity scores are
essential in order to balance covariates between internal and external
trial participants.

**Usage:**

``` r
calc_prop_scr(internal_df, external_df, id_col, model)
```

- `internal_df`: Internal dataset.
- `external_df`: External dataset.
- `id_col`: Name of the column in both datasets used to identify each
  subject.
- `model`: Model used to calculate propensity scores.

### `calc_power_prior()`

Calculate power priors based on a distribution, hyperparameters, a
weighted object and a response variable. Power priors provide a Bayesian
framework for incorporating external information into the analysis.

**Usage:**

``` r
calc_power_prior(prior, weighted_obj, response)
```

- `prior`: a {distributional} object that is the prior of the external
  data.
- `weighted_obj`: A weighted object created by calling
  `calc_prop_scr()`.
- `response`: Name of the response variable.
- `...`: Additional parameters needed depending on the type of power
  prior used.

## Examples

Here are some examples demonstrating the usage of the package:

``` r
library(beastt)
library(ggdist)
library(ggplot2)

set.seed(1234)

# Create internal and external datasets
internal_df <- data.frame(id_col = 1:20, cov1 = rnorm(10, 2), cov2 = rnorm(100, 20), 
                          y = rnorm(20, mean = 5, sd = 3)
external_df <- data.frame(id_col = 21:40, cov1 = rnorm(10, 2), cov2 = rnorm(100, 18), 
                          y = rnorm(20, mean = 8, sd = 4)

# Example for propensity score calculation
model <- as.formula("~cov1 + cov2")
psscores <- calc_prop_scr(internal_df = internal_df, external_df = external_df, 
                           id_col = id_col, model = model)

# Example for power prior calculation using a Normal prior for external data
powerprior <- calc_power_prior(prior = dist_normal(15, 100), weighted_obj = psscores, 
                               response = y)

# Plot power prior and interval
ggplot(tibble(), aes(xdist = powerprior, y = 1)) +
  stat_slabinterval()
```

## Contributing

Feel free to contribute to the “beastt” package by reporting issues or
submitting pull requests on the GitHub repository.

## License

This package is released under the [Apache License](%3E=%202).
