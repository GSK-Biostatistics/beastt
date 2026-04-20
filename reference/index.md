# Package index

## Propensity Scores

- [`calc_prop_scr()`](https://gsk-biostatistics.github.io/beastt/reference/calc_prop_scr.md)
  : Create a Propensity Score Object

- [`trim_ps()`](https://gsk-biostatistics.github.io/beastt/reference/trim_ps.md)
  :

  Trim a `prop_scr` object

- [`rescale_ps()`](https://gsk-biostatistics.github.io/beastt/reference/rescale_ps.md)
  :

  Rescale a `prop_scr` object

- [`tidy(`*`<prop_scr>`*`)`](https://gsk-biostatistics.github.io/beastt/reference/tidy.prop_scr.md)
  : Tidy a(n) prop_scr object

- [`is_prop_scr()`](https://gsk-biostatistics.github.io/beastt/reference/is_prop_scr.md)
  : Test If Propensity Score Object

## Power Prior

- [`calc_power_prior_beta()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_beta.md)
  : Calculate Power Prior Beta
- [`calc_power_prior_norm()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_norm.md)
  : Calculate Power Prior Normal
- [`calc_power_prior_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_power_prior_weibull.md)
  : Calculate Power Prior Weibull

## Robustify

- [`robustify_mvnorm()`](https://gsk-biostatistics.github.io/beastt/reference/robustify_mvnorm.md)
  : Robustify Multivariate Normal Distributions
- [`robustify_norm()`](https://gsk-biostatistics.github.io/beastt/reference/robustify_norm.md)
  : Robustify Normal Distributions

## Posterior

- [`calc_post_beta()`](https://gsk-biostatistics.github.io/beastt/reference/calc_post_beta.md)
  : Calculate Posterior Beta
- [`calc_post_norm()`](https://gsk-biostatistics.github.io/beastt/reference/calc_post_norm.md)
  : Calculate Posterior Normal
- [`calc_post_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_post_weibull.md)
  : Calculate Posterior Weibull

## Simulation

- [`bootstrap_cov()`](https://gsk-biostatistics.github.io/beastt/reference/bootstrap_cov.md)
  : Bootstrap Covariate Data
- [`calc_cond_binary()`](https://gsk-biostatistics.github.io/beastt/reference/calc_cond_binary.md)
  : Calculate Conditional Drift and Treatment Effect for Binary Outcome
  Models
- [`inv_logit()`](https://gsk-biostatistics.github.io/beastt/reference/inv_logit.md)
  : Inverse Logit Function
- [`sim_accrual()`](https://gsk-biostatistics.github.io/beastt/reference/sim_accrual.md)
  : Simulate Participant Accrual Times
- [`sim_weib_ph()`](https://gsk-biostatistics.github.io/beastt/reference/sim_weib_ph.md)
  : Simulate Event Times for Each Participant from a Weibull
  Proportional Hazards Regression Model
- [`sim_pw_const_haz()`](https://gsk-biostatistics.github.io/beastt/reference/sim_pw_const_haz.md)
  : Simulate Event Times for Each Individual from a Piecewise Constant
  Hazard Model
- [`calc_cond_weibull()`](https://gsk-biostatistics.github.io/beastt/reference/calc_cond_weibull.md)
  : Calculate Conditional Drift and Treatment Effect for Time-to-Event
  Outcome Models
- [`calc_study_duration()`](https://gsk-biostatistics.github.io/beastt/reference/calc_study_duration.md)
  : Calculate the Analysis Time Based on a Target Number of Events
  and/or Target Follow-up Time

## Visualisation

- [`prop_scr_cloud()`](https://gsk-biostatistics.github.io/beastt/reference/prop_scr_cloud.md)
  : Propensity Score Cloud Plot
- [`prop_scr_dens()`](https://gsk-biostatistics.github.io/beastt/reference/prop_scr_dens.md)
  : Density of the Propensity Score Object
- [`prop_scr_hist()`](https://gsk-biostatistics.github.io/beastt/reference/prop_scr_hist.md)
  : Histogram of the Propensity Score Object
- [`prop_scr_love()`](https://gsk-biostatistics.github.io/beastt/reference/prop_scr_love.md)
  : Love Plot of the Absolute Standardized Mean Differences
- [`plot_dist()`](https://gsk-biostatistics.github.io/beastt/reference/plot_dist.md)
  : Plot Distribution
- [`sweet_spot_plot()`](https://gsk-biostatistics.github.io/beastt/reference/sweet_spot_plot.md)
  : Create Sweet Spot Plots for Multiple Simulation Scenarios

## Utility

- [`mix_means()`](https://gsk-biostatistics.github.io/beastt/reference/mix_means.md)
  : Extract Means of Mixture Components
- [`mix_sigmas()`](https://gsk-biostatistics.github.io/beastt/reference/mix_sigmas.md)
  : Extract Standard Deviations of Mixture Components
- [`avg_dist()`](https://gsk-biostatistics.github.io/beastt/reference/avg_dist.md)
  : Calculate Average Distribution from Multiple Distributional Objects
- [`approx_mvn_at_time()`](https://gsk-biostatistics.github.io/beastt/reference/approx_mvn_at_time.md)
  : Approximate Multivariate Normal Distribution as Beta at a Specific
  Time
- [`beastt-package`](https://gsk-biostatistics.github.io/beastt/reference/beastt-package.md)
  [`beastt`](https://gsk-biostatistics.github.io/beastt/reference/beastt-package.md)
  : The 'beastt' package.

## Data

- [`ex_binary_df`](https://gsk-biostatistics.github.io/beastt/reference/ex_binary_df.md)
  : External Binary Control Data for Propensity Score Balancing
- [`ex_norm_df`](https://gsk-biostatistics.github.io/beastt/reference/ex_norm_df.md)
  : External Normal Control Data for Propensity Score Balancing
- [`ex_tte_df`](https://gsk-biostatistics.github.io/beastt/reference/ex_tte_df.md)
  : External Time-to-Event Control Data for Propensity Score Balancing
- [`int_binary_df`](https://gsk-biostatistics.github.io/beastt/reference/int_binary_df.md)
  : Internal Binary Data for Propensity Score Balancing
- [`int_norm_df`](https://gsk-biostatistics.github.io/beastt/reference/int_norm_df.md)
  : Internal Normal Data for Propensity Score Balancing
- [`int_tte_df`](https://gsk-biostatistics.github.io/beastt/reference/int_tte_df.md)
  : Internal Time-to-Event Control Data for Propensity Score Balancing
- [`binary_sim_df`](https://gsk-biostatistics.github.io/beastt/reference/binary_sim_df.md)
  : Binary Simulation Data
- [`tte_sim_df`](https://gsk-biostatistics.github.io/beastt/reference/tte_sim_df.md)
  : Time-to-Event Simulation Data
