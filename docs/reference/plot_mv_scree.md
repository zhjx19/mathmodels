# Scree Plot

Plots the proportion and cumulative proportion of variance explained by
each component (PCA) or factor (EFA).

## Usage

``` r
plot_mv_scree(x)
```

## Arguments

- x:

  An `mv_pca` or `mv_efa` result.

## Value

A `ggplot` object.

## Examples

``` r
r = mv_pca(mtcars, cyl:carb)
plot_mv_scree(r)
```
