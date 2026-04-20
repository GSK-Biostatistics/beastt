# Trim a `prop_scr` object

Trim a `prop_scr` object

## Usage

``` r
trim_ps(x, low = NULL, high = NULL, quantile = FALSE)
```

## Arguments

- x:

  A `prop_scr` object

- low:

  Low cut-off such that all participants with propensity scores less
  than this value (or quantile if `quantile = TRUE`) are removed. If
  left `NULL` no lower bound will be used

- high:

  High cut-off such that all participants with propensity scores greater
  than this value (or quantile if `quantile = TRUE`) are removed. If
  left `NULL` no upper bound will be used

- quantile:

  True/False value to determine if the cut-off values are based directly
  on the propensity scores (false) or their quantiles (true). By default
  this is false.

## Value

a `prop_scr` object with a trimmed propensity score distribution

## Details

This function uses R's default method of quantile calculation (type 7)

## Examples

``` r
library(dplyr)
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                       external_df = ex_binary_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
trim_ps(ps_obj, low = 0.3, high = 0.7)
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> • cov1 + cov2 + cov3 + cov4
#> 
#> ── Propensity Scores and Weights ───────────────────────────────────────────────
#> • Effective sample size of the external arm: 60
#> # A tibble: 84 × 4
#>    subjid Internal `Propensity Score` `Inverse Probability Weight`
#>     <int> <lgl>                 <dbl>                        <dbl>
#>  1      1 FALSE                 0.333                        0.500
#>  2      3 FALSE                 0.539                        1.17 
#>  3      4 FALSE                 0.546                        1.20 
#>  4      5 FALSE                 0.344                        0.524
#>  5      6 FALSE                 0.393                        0.646
#>  6      7 FALSE                 0.390                        0.639
#>  7      8 FALSE                 0.340                        0.515
#>  8     12 FALSE                 0.356                        0.553
#>  9     15 FALSE                 0.344                        0.524
#> 10     17 FALSE                 0.440                        0.785
#> # ℹ 74 more rows
#> 
#> ── Absolute Standardized Mean Difference ───────────────────────────────────────
#> # A tibble: 4 × 3
#>   covariate diff_unadj diff_adj
#>   <chr>          <dbl>    <dbl>
#> 1 cov1          0.265    0.434 
#> 2 cov2          0.0367   0.0163
#> 3 cov3          0.0123   0.0950
#> 4 cov4          0.272    0.372 
```
