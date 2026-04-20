# Inverse Logit Function

Inverse Logit Function

## Usage

``` r
inv_logit(x)
```

## Arguments

- x:

  Real number(s) to take the inverse logit of

## Value

Vector of probabilities between 0 and 1

## Details

This function is a short hand to \\\exp(x)/(1 + \exp(x))\\.

## Examples

``` r
inv_logit(0.5)
#> [1] 0.6224593
```
