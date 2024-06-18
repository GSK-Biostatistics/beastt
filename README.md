
# beastt

## Bayesian Evaluation, Analysis, and Simulation Software Tools for Trials (beastt)

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/beastt)](https://CRAN.R-project.org/package=beastt)
[![R-CMD-check](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

Welcome to the “beastt” package! This R package is designed to assist
users in performing Bayesian dynamic borrowing with covariate adjustment
via inverse probability weighting for simulations and data analyses in
clinical trials.

## Installation

You can install the development version of {beastt} from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GSK-Biostatistics/beastt")
```

## Usage

At the moment {beastt} covers cases when borrowing from external control
data with either normal or binary endpoints. For more information, see
the vignettes. Future updates are expected to include cases with
survival endpoints and repeated measure.

## Contributing

Feel free to contribute to the {beastt} package by reporting issues or
submitting pull requests on the GitHub repository.

## License

This package is released under the Apache License.
