
# beastt

## Bayesian Evaluation, Analysis, and Simulation Software Tools for Trials (beastt)

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/beastt)](https://CRAN.R-project.org/package=beastt)
[![R-CMD-check](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GSK-Biostatistics/beastt/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/GSK-Biostatistics/beastt/graph/badge.svg?token=UVH8SF4OXT)](https://app.codecov.io/gh/GSK-Biostatistics/beastt)

<!-- badges: end -->

## Overview

Welcome to the {[beastt](https://gsk-biostatistics.github.io/beastt/)}
package! This R package is designed to assist users in performing
Bayesian dynamic borrowing with covariate adjustment via inverse
probability weighted robust mixture priors for simulations and data
analyses in clinical trials. For the sake of this package, we use the
term IPW BMB to refer to this inverse probability weighted robust
mixture methodology.

## What is IPW BDB?

Inverse Probability Weighted Bayesian Dynamic Borrowing (IPW BDB) is a
statistical approach designed to enhance the estimation of marginal
(i.e., population-averaged or unconditional) treatment effects in
clinical trials. This method employs inverse probability weighted robust
mixture priors to adjust for covariate differences between a new
internal study and external (i.e., historical) control data.

## Why use IPW BDB?

By using propensity score-based inverse probability weighting, IPW BDB
effectively balances prognostic variables between trial participants and
historical controls, improving inference accuracy and reducing biases
due to differences in covariate distributions. This technique increases
the statistical power and reduces potential biases in estimating average
treatment effects, which are critical for health policy decisions and
drug approval processes.

IPW BDB has two mechanisms by which it can account for drift from
different sources: 1. The use of a robust mixture prior alleviates
prior-data conflict by dynamically down weighting external data when
there is a significant level of drift between studies. 2. Inverse
probability weighting can account for explainable causes of drift by
balancing covariate distributions between external and internal control
participants.

*Augmenting the standard robust mixture prior (RMP) approach to
incorporate IPWs does not add substantial computational burden
associated with other Bayesian approaches*; e.g., in cases where
conjugate priors exist for the standard RMP approach, they will still
exist for the IPW BDB approach.

## When can IPW BDB be used?

IPW BDB should be considered in clinical trial settings where individual
level external control data is available and you want to integrate this
data with your current trial. This method is particularly useful when
there are differences in distributions of key prognostic factors between
the current study population and the external controls, which could
otherwise introduce bias. It is especially relevant in oncology and rare
disease trials, where using external data can help overcome challenges
such as slow patient enrollment due to reluctance to join control
groups. IPW BDB is well-suited for contexts where Bayesian dynamic
borrowing is already applicable but could benefit from additional
adjustments for confounding.

## Installation

You can install the development version of {beastt} from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GSK-Biostatistics/beastt")
```

## Usage

At the moment {beastt} covers borrowing from external control data for
normal, binary, and time to event endpoints. For more information, see
the vignettes.

## Contributing

Feel free to contribute to the {beastt} package by reporting issues or
submitting pull requests on the GitHub repository.

## License

This package is released under GLP-3.
