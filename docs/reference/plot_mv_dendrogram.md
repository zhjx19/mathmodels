# Dendrogram Plot

Draws the dendrogram of an `mv_hclust` result with ggplot2, with a
dashed reference line at the tree height where the `k` clusters are cut.

## Usage

``` r
plot_mv_dendrogram(x)
```

## Arguments

- x:

  An `mv_hclust` result.

## Value

A `ggplot` object.

## Examples

``` r
r = mv_hclust(mtcars, cyl:carb, k = 3)
plot_mv_dendrogram(r)
```
