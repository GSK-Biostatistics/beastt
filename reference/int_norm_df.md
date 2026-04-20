# Internal Normal Data for Propensity Score Balancing

This is a simulated dataset used to illustrate Bayesian dynamic
borrowing in the case when borrowing from an external control arm with a
normal endpoint, where the baseline covariate distributions of the
internal and external data are balanced via inverse probability
weighting.

## Usage

``` r
int_norm_df
```

## Format

### `int_norm_df`

A data frame with 120 rows and 7 columns:

- subjid:

  Unique subject ID

- cov1:

  Covariate 1, which is normally distributed around 55 with a SD of 8

- cov2:

  Covariate 2, which is binary (0 vs. 1) with about 30% of participants
  having level 1

- cov3:

  Covariate 3, which is binary (0 vs. 1) with about 50% of participants
  having level 1

- cov4:

  Covariate 4, which is binary (0 vs. 1) with about 30% of participants
  having level 1

- trt:

  Treatment indicator, where 0 = control and 1 = active treatment

- y:

  Response, which is normally distributed with a SD of 0.15
