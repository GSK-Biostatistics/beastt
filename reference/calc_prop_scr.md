# Create a Propensity Score Object

Calculate the propensity scores and ATT inverse probability weights for
participants from internal and external datasets. Only the relevant
treatment arms from each dataset should be read in (e.g., only the
control arm from each dataset if creating a hybrid control arm).

## Usage

``` r
calc_prop_scr(internal_df, external_df, id_col, model, ...)
```

## Arguments

- internal_df:

  Internal dataset with one row per subject and all the variables needed
  to run the model

- external_df:

  External dataset with one row per subject and all the variables needed
  to run the model

- id_col:

  Name of the column in both datasets used to identify each subject. It
  must be the same across datasets

- model:

  Model used to calculate propensity scores

- ...:

  Optional arguments

## Value

`prop_scr_obj` object, with the internal and the external data and the
propensity score and inverse probability weight calculated for each
subject.

## Details

For the subset of participants in both the external and internal studies
for which we want to balance the covariate distributions (e.g., external
control and internal control participants if constructing a hybrid
control arm), we define a study-inclusion propensity score for each
participant as

\$\$e(x_i) = P(S_i = 1 \mid x_i),\$\$

where \\x_i\\ denotes a vector of baseline covariates for the \\i\\th
participant and \\S_i\\ denotes the indicator that the participant is
enrolled in the internal trial (\\S_i = 1\\ if internal, \\S_i = 0\\ if
external). The estimated propensity score \\\hat{e}(x_i)\\ is obtained
using logistic regression.

An ATT inverse probability weight is calculated for each individual as

\$\$\hat{a}\_{0i} = \frac{\hat{e}(x_i)}{\hat{P}(S_i = s_i \| x_i)} =
s_i + (1 - s_i ) \frac{\hat{e}(x_i)}{1 - \hat{e}(x_i)}.\$\$

In a weighted estimator, data from participants in the external study
are given a weight of \\\hat{e}(x_i)⁄(1 - \hat{e}(x_i))\\ whereas data
from participants in the internal trial are given a weight of 1.

## Examples

``` r
# This can be used for both continuous and binary data
library(dplyr)
# Continuous
calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                       external_df = ex_norm_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> • cov1 + cov2 + cov3 + cov4
#> 
#> ── Propensity Scores and Weights ───────────────────────────────────────────────
#> • Effective sample size of the external arm: 61
#> # A tibble: 150 × 4
#>    subjid Internal `Propensity Score` `Inverse Probability Weight`
#>     <int> <lgl>                 <dbl>                        <dbl>
#>  1      1 FALSE                 0.181                        0.221
#>  2      2 FALSE                 0.351                        0.542
#>  3      3 FALSE                 0.456                        0.840
#>  4      4 FALSE                 0.115                        0.130
#>  5      5 FALSE                 0.364                        0.572
#>  6      6 FALSE                 0.377                        0.604
#>  7      7 FALSE                 0.131                        0.150
#>  8      8 FALSE                 0.259                        0.349
#>  9      9 FALSE                 0.235                        0.308
#> 10     10 FALSE                 0.207                        0.261
#> # ℹ 140 more rows
#> 
#> ── Absolute Standardized Mean Difference ───────────────────────────────────────
#> # A tibble: 4 × 3
#>   covariate diff_unadj diff_adj
#>   <chr>          <dbl>    <dbl>
#> 1 cov1          0.506   0.104  
#> 2 cov2          0.0548  0.104  
#> 3 cov3          0.0667  0.00115
#> 4 cov4          0.275   0.00266
# Binary
calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                       external_df = ex_binary_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> • cov1 + cov2 + cov3 + cov4
#> 
#> ── Propensity Scores and Weights ───────────────────────────────────────────────
#> • Effective sample size of the external arm: 81
#> # A tibble: 150 × 4
#>    subjid Internal `Propensity Score` `Inverse Probability Weight`
#>     <int> <lgl>                 <dbl>                        <dbl>
#>  1      1 FALSE                 0.333                        0.500
#>  2      2 FALSE                 0.288                        0.405
#>  3      3 FALSE                 0.539                        1.17 
#>  4      4 FALSE                 0.546                        1.20 
#>  5      5 FALSE                 0.344                        0.524
#>  6      6 FALSE                 0.393                        0.646
#>  7      7 FALSE                 0.390                        0.639
#>  8      8 FALSE                 0.340                        0.515
#>  9      9 FALSE                 0.227                        0.294
#> 10     10 FALSE                 0.280                        0.389
#> # ℹ 140 more rows
#> 
#> ── Absolute Standardized Mean Difference ───────────────────────────────────────
#> # A tibble: 4 × 3
#>   covariate diff_unadj diff_adj
#>   <chr>          <dbl>    <dbl>
#> 1 cov1          0.339  0.0461  
#> 2 cov2          0.0450 0.0204  
#> 3 cov3          0.160  0.000791
#> 4 cov4          0.308  0.00857 
```
