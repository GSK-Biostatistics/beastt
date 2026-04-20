# Simulate Event Times for Each Individual from a Piecewise Constant Hazard Model

Simulate Event Times for Each Individual from a Piecewise Constant
Hazard Model

## Usage

``` r
sim_pw_const_haz(n, hazard_periods = NULL, hazard_values)
```

## Arguments

- n:

  Number of individuals

- hazard_periods:

  Vector of break points between time periods with separate constant
  hazards, e.g., c(6,8) defines \[0,6), \[6,8), \[8, infinity). Leave as
  NULL if defining only one hazard period.

- hazard_values:

  Vector of constant hazard values associated with the time intervals

## Value

Vector of simulated times from the time-to-event distribution

## Details

Simulate the event or censoring times for each participant using a
piecewise constant hazard model where each time period is defined to
have a different constant hazard. Here, `hazard_periods` defines the
right time points for each time period (with the exception of the last
time period which extends to infinity), and `hazard_values` defines the
constant hazards for each time period.

## Examples

``` r
tte_dat <- sim_pw_const_haz(n = 100000, hazard_periods = c(6, 8), hazard_values = c(0.1, 0.1, 0.1))
hist(tte_dat, breaks = 100, main = "Event Time Distribution", xlab = "Event Time")
```
