# External Normal Control Data for Propensity Score Balancing

This is a simulated dataset used to illustrate Bayesian dynamic
borrowing in the case when borrowing from an external control arm with a
normal endpoint, where the baseline covariate distributions of the
internal and external data are balanced via inverse probability
weighting.

## Usage

``` r
ex_norm_df
```

## Format

### `ex_norm_df`

A data frame with 150 rows and 6 columns:

- subjid:

  Unique subject ID

- cov1:

  Covariate 1, which is normally distributed around 50 with a SD of 10

- cov2:

  Covariate 2, which is binary (0 vs. 1) with about 20% of participants
  having level 1

- cov3:

  Covariate 3, which is binary (0 vs. 1) with about 60% of participants
  having level 1

- cov4:

  Covariate 4, which is binary (0 vs. 1) with about 30% of participants
  having level 1

- y:

  Response, which is normally distributed with a SD of 0.15
