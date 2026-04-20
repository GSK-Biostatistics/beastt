# Plot Distribution

Plot Distribution

## Usage

``` r
plot_dist(...)
```

## Arguments

- ...:

  Distributional object(s) to plot. When passing multiple objects naming
  them will change the labels in the plot, else they will use the
  distributional format

## Value

ggplot object that is the density of the provided distribution

## Examples

``` r
library(distributional)
plot_dist(dist_normal(0, 1))

plot_dist(dist_multivariate_normal(mu = list(c(1, 2)), sigma = list(matrix(c(4, 2, 2, 3), ncol=2))))

#Plotting Multiple
plot_dist(dist_normal(0, 1), dist_normal(10, 5))

plot_dist('Prior' = dist_normal(0, 1), 'Posterior' = dist_normal(10, 5))
```
