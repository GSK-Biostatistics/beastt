# Time-to-Event Simulation Data

This is an example of output from a simulation study that investigates
the operating characteristics of inverse probability weighted Bayesian
dynamic borrowing for the case with a time-to-event outcome. This output
was generated based on the time-to-event simulation template. For this
simulation study, only the degree of covariate imbalance (as indicated
by `population`) and the marginal treatment effect were varied.

## Usage

``` r
tte_sim_df
```

## Format

### `tte_sim_df` A data frame with 18 rows and 7 columns:

- population:

  Populations defined by different covariate imbalances

- marg_trt_eff:

  Marginal treatment effect

- true_control_surv_prop:

  True control survival probability at time t=12 months on the marginal
  scale

- reject_H0_yes:

  Probability of rejecting the null hypothesis in the case with
  borrowing

- no_borrowing_reject_H0_yes:

  Probability of rejecting the null hypothesis without borrowing

- pwr_prior:

  Vector of IPW power priors as distributional objects

- mix_prior:

  Vector of mixture priors (i.e., the robustified IPW power priors) as
  distributional objects
