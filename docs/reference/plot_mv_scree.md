# Scree Plot

Plots the proportion or cumulative proportion of variance explained by
each component (PCA) or factor (EFA).

## Usage

``` r
plot_mv_scree(x, type = c("prop", "cum"))
```

## Arguments

- x:

  An `mv_pca` or `mv_efa` result.

- type:

  `"prop"` (default) draws only the per-component proportion of
  variance; `"cum"` draws only the cumulative proportion.

## Value

A `ggplot` object.

## Examples

``` r
r = mv_pca(mtcars, cyl:carb)
plot_mv_scree(r)
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?

plot_mv_scree(r, type = "cum")
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
```
