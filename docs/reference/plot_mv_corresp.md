# Correspondence Analysis Map

Draws the symmetric two-dimensional map of a `mv_corresp` result: row
and column coordinates with labels.

## Usage

``` r
plot_mv_corresp(x)
```

## Arguments

- x:

  An `mv_corresp` result with at least two dimensions.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(42)
df = data.frame(
  breed = rep(c("A", "B", "C"), each = 30),
  grade = sample(c("low", "mid", "high"), 90, TRUE)
)
r = mv_corresp(df, breed, grade)
plot_mv_corresp(r)
```
