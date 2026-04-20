# Rescale a `prop_scr` object

Rescale a `prop_scr` object

## Usage

``` r
rescale_ps(x, n = NULL, scale_factor = NULL)
```

## Arguments

- x:

  a `prop_scr` obj

- n:

  Desired sample size that the external data should effectively
  contribute to the analysis of the internal trial data. This will be
  used to scale the external weights if `scale_factor` is not specified

- scale_factor:

  Value to multiple all weights by. This will be used to scale the
  external weights if `n` is not specified

## Value

a `prop_scr` object with rescaled weights

## Examples

``` r
library(dplyr)
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                       external_df = ex_binary_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
# weights in a propensity score object can be rescaled to achieve a desired
# effective sample size (i.e., sum of weights)
rescale_ps(ps_obj, n = 75)
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> • cov1 + cov2 + cov3 + cov4
#> 
#> ── Propensity Scores and Weights ───────────────────────────────────────────────
#> • Effective sample size of the external arm: 75
#> # A tibble: 150 × 4
#>    subjid Internal `Propensity Score` `Inverse Probability Weight`
#>     <int> <lgl>                 <dbl>                        <dbl>
#>  1      1 FALSE                 0.333                        0.465
#>  2      2 FALSE                 0.288                        0.377
#>  3      3 FALSE                 0.539                        1.09 
#>  4      4 FALSE                 0.546                        1.12 
#>  5      5 FALSE                 0.344                        0.487
#>  6      6 FALSE                 0.393                        0.601
#>  7      7 FALSE                 0.390                        0.594
#>  8      8 FALSE                 0.340                        0.479
#>  9      9 FALSE                 0.227                        0.273
#> 10     10 FALSE                 0.280                        0.362
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

# Or by a predetermined factor
rescale_ps(ps_obj, scale_factor = 1.5)
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> • cov1 + cov2 + cov3 + cov4
#> 
#> ── Propensity Scores and Weights ───────────────────────────────────────────────
#> • Effective sample size of the external arm: 121
#> # A tibble: 150 × 4
#>    subjid Internal `Propensity Score` `Inverse Probability Weight`
#>     <int> <lgl>                 <dbl>                        <dbl>
#>  1      1 FALSE                 0.333                        0.751
#>  2      2 FALSE                 0.288                        0.608
#>  3      3 FALSE                 0.539                        1.75 
#>  4      4 FALSE                 0.546                        1.81 
#>  5      5 FALSE                 0.344                        0.785
#>  6      6 FALSE                 0.393                        0.970
#>  7      7 FALSE                 0.390                        0.958
#>  8      8 FALSE                 0.340                        0.773
#>  9      9 FALSE                 0.227                        0.441
#> 10     10 FALSE                 0.280                        0.584
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
