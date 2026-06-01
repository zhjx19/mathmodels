# Compute fuzzy membership vector and return corresponding membership functions.

`compute_mf` transforms a single indicator value into a fuzzy membership
vector, where each element represents the degree of membership to a
specific evaluation level. `compute_mf_funs` returns the list of
membership functions for visualization purposes.

## Usage

``` r
compute_mf_funs(thresholds)

compute_mf(x, thresholds)
```

## Arguments

- thresholds:

  A numeric vector containing at least two threshold values that define
  the boundaries between evaluation levels.

- x:

  A numeric scalar representing the value of an indicator.

## Value

A list with two elements:

- mv:

  A numeric vector, membership degrees to each level.

- mfs:

  A list of functions, one per level, for plotting membership functions.

## Examples

``` r
# Example: SO2 concentration = 0.07, thresholds = c(0.05, 0.15, 0.25, 0.5)
th = c(0.05, 0.15, 0.25, 0.5)
compute_mf(0.07, th)
#> [1] 0.8 0.2 0.0 0.0

if (FALSE) { # \dontrun{
mfs = compute_mf_funs(th)
plots = lapply(mfs, \(x) plot_mf(x, xlim = c(0, 0.6)))
gridExtra::grid.arrange(grobs = plots, nrow = 2)
} # }

```
