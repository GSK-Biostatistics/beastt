# Test If Propensity Score Object

Test If Propensity Score Object

## Usage

``` r
is_prop_scr(x)
```

## Arguments

- x:

  Object to test

## Value

Boolean

## Examples

``` r
library(dplyr)
x <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                       external_df = ex_norm_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
is_prop_scr(x)
#> [1] TRUE
```
