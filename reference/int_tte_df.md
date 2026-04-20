# Internal Time-to-Event Control Data for Propensity Score Balancing

This is a simulated dataset used to illustrate Bayesian dynamic
borrowing in the case when borrowing from an external control arm with a
time-to-event endpoint, where the baseline covariate distributions of
the internal and external data are balanced via inverse probability
weighting.

## Usage

``` r
int_tte_df
```

## Format

### `int_tte_df`

A data frame with 160 rows and 10 columns:

- subjid:

  Unique subject ID

- y:

  Response (observed time at which the participant either had an event
  or was censored)

- enr_time:

  Enrollment time

- total_time:

  Time from study start

- event:

  Event indicator (1: event; 0: censored)

- trt:

  Treatment indicator, where 0 = control and 1 = active treatment

- cov1:

  Covariate 1, which is normally distributed around 62 with a SD of 8

- cov2:

  Covariate 2, which is binary (0 vs. 1) with about 40% of participants
  having level 1

- cov3:

  Covariate 3, which is binary (0 vs. 1) with about 40% of participants
  having level 1

- cov4:

  Covariate 4, which is binary (0 vs. 1) with about 60% of participants
  having level 1
