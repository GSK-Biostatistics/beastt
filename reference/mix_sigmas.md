# Extract Standard Deviations of Mixture Components

Extract Standard Deviations of Mixture Components

## Usage

``` r
mix_sigmas(x)
```

## Arguments

- x:

  A mixture distributional object

## Value

numeric or list object

## Details

If a distributional object that is a mixture of two or more normal
distributions is read in, the function will return a numeric object with
the standard deviations of each normal component. If the distributional
object is a mixture of two or more multivariate normal distributions,
the function will return a list with the covariance matrices of each
multivariate normal component.

## Examples

``` r
library(distributional)
mix_norm <- dist_mixture(comp1 = dist_normal(1, 10),
                         comp2 = dist_normal(1.5, 12),
                         weights = c(.5, .5))
mix_sigmas(mix_norm)
#> comp1 comp2 
#>    10    12 
```
