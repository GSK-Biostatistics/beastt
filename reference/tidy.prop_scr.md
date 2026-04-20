# Tidy a(n) prop_scr object

Tidy a(n) prop_scr object

## Usage

``` r
# S3 method for class 'prop_scr'
tidy(x, ...)
```

## Arguments

- x:

  a `prop_scr` obj

- ...:

  Unused, included for generic consistency only.

## Value

A tidy
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
summarizing the results of the propensity score weighting. The tibble
will have the id column of the external data, an `internal` column to
indicate all the data is external, a `ps` column with the propensity
scores and a `weight` column with the inverse probability weights

## Examples

``` r
library(dplyr)
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                       external_df = ex_binary_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
tidy(ps_obj)
#> # A tibble: 150 × 4
#>    subjid internal    ps weight
#>     <int> <lgl>    <dbl>  <dbl>
#>  1      1 FALSE    0.333  0.500
#>  2      2 FALSE    0.288  0.405
#>  3      3 FALSE    0.539  1.17 
#>  4      4 FALSE    0.546  1.20 
#>  5      5 FALSE    0.344  0.524
#>  6      6 FALSE    0.393  0.646
#>  7      7 FALSE    0.390  0.639
#>  8      8 FALSE    0.340  0.515
#>  9      9 FALSE    0.227  0.294
#> 10     10 FALSE    0.280  0.389
#> # ℹ 140 more rows
```
