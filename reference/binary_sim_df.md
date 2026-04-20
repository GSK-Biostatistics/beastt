# Binary Simulation Data

This is an example of output from a simulation study that investigates
the operating characteristics of inverse probability weighted Bayesian
dynamic borrowing for the case with a binary outcome. This output was
generated based on the binary simulation template. For this simulation
study, only the degree of covariate imbalance (as indicated by
`population`) and the marginal treatment effect were varied.

## Usage

``` r
binary_sim_df
```

## Format

### `binary_sim_df` A data frame with 255 rows and 6 columns:

- population:

  Populations defined by different covariate imbalances

- marg_trt_eff:

  Marginal treatment effect

- true_control_RR:

  True control response rate on the marginal scale

- reject_H0_yes:

  Probability of rejecting the null hypothesis in the case with
  borrowing

- no_borrowing_reject_H0_yes:

  Probability of rejecting the null hypothesis without borrowing

- pwr_prior:

  Vector of power priors (or some other informative prior distribution
  for the control marginal parameter of interest based on the external
  data) as distributional objects
