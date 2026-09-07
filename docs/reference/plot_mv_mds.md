# MDS Configuration Plot

Plots the first two dimensions of an `mv_mds` configuration.

## Usage

``` r
plot_mv_mds(x)
```

## Arguments

- x:

  An `mv_mds` result with `k >= 2`.

## Value

A `ggplot` object.

## Examples

``` r
r = mv_mds(mtcars, cyl:carb, k = 2)
plot_mv_mds(r)
```
