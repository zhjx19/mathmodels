# Scores Plot

Plots the first two score dimensions of an `mv_pca` result (PC1 vs PC2)
or an LDA `mv_discrim` result (LD1 vs LD2), optionally colored by group.

## Usage

``` r
plot_mv_scores(x, group = NULL)
```

## Arguments

- x:

  An `mv_pca` result, or an `mv_discrim` result from LDA (QDA has no
  scores); LDA with only 2 groups yields a single discriminant function,
  which is not enough for a 2-D scores plot.

- group:

  Optional vector with one value per row of the scores, used for
  coloring.

## Value

A `ggplot` object.

## Examples

``` r
r = mv_pca(mtcars, cyl:carb)
plot_mv_scores(r)
```
