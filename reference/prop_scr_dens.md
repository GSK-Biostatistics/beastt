# Density of the Propensity Score Object

Plot overlapping density curves of the propensity scores for both the
internal and external participants, or plot external IPWs.

## Usage

``` r
prop_scr_dens(
  x,
  variable = c("propensity score", "ps", "inverse probability weight", "ipw"),
  ...
)
```

## Arguments

- x:

  Propensity score object

- variable:

  Variable to plot. It must be either a propensity score ("ps" or
  "propensity score") or inverse probability weight ("ipw" or "inverse
  probability weight")

- ...:

  Optional arguments for `geom_density`

## Value

ggplot object

## Examples

``` r
library(dplyr)
ps_obj <- calc_prop_scr(internal_df = filter(int_norm_df, trt == 0),
                       external_df = ex_norm_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
# Plotting the Propensity Scores
prop_scr_dens(ps_obj)

# Or plotting the inverse probability weights
prop_scr_dens(ps_obj, variable = "ipw")

```
