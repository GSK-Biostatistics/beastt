# Simulate Participant Accrual Times

Simulate Participant Accrual Times

## Usage

``` r
sim_accrual(n, accrual_periods, accrual_props)
```

## Arguments

- n:

  Number of participants

- accrual_periods:

  Vector of right endpoints defining the time periods of accrual, e.g.,
  c(6,8) defines 0\<=x\<6, 6\<=x\<8.

- accrual_props:

  Vector indicating the proportion of participants that are expected to
  be enrolled during each of the defined accrual periods. Should sum to
  1, otherwise these proportions will be normalized.

## Value

Vector of accrual times corresponding to when each participant enters
the study

## Details

Simulate the accrual times for each participant, where `accrual_periods`
defines the right time points for each accrual period (with the last
element corresponding to the end of accrual), and `accrual_props`
defines the proportion of study participants who are enrolled during
each of these periods. The simulated accrual times for participants
within a given accrual period are assumed to be uniformly distributed.

## Examples

``` r
at <- sim_accrual(n = 100000, accrual_periods = c(6, 8), accrual_props = c(.5, .5))
hist(at, breaks = 100, main = "Histogram of Enrollment Times", xlab = "Enrollment Time")
```
