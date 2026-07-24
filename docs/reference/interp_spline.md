# Cubic Spline Interpolation

Performs cubic spline interpolation. When `spar` is `NULL`, uses
[`stats::spline()`](https://rdrr.io/r/stats/splinefun.html) with the FMM
method. Otherwise uses
[`stats::smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html)
for smoothing.

## Usage

``` r
interp_spline(x, y, xout, spar = NULL)
```

## Arguments

- x:

  Numeric vector of observed x-values (no duplicates, length \>= 4).

- y:

  Numeric vector of observed y-values.

- xout:

  Numeric vector of x-values where interpolation is desired.

- spar:

  Smoothing parameter passed to
  [`smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html).
  `NULL` for standard (non-smoothing) spline.

## Value

A tibble with columns `x` and `y`. The attribute `method` is set to
`"spline"` or `"smooth_spline"` depending on the variant used.

## Examples

``` r
x = c(1, 2, 3, 4, 5)
y = c(2, 4, 1, 5, 3)
interp_spline(x, y, xout = c(1.5, 2.5, 3.5))
#> # A tibble: 3 × 2
#>       x     y
#>   <dbl> <dbl>
#> 1   1.5  4.55
#> 2   2.5  2.03
#> 3   3.5  2.59
interp_spline(x, y, xout = c(1.5, 2.5, 3.5), spar = 0.8)
#> # A tibble: 3 × 2
#>       x     y
#>   <dbl> <dbl>
#> 1   1.5  2.55
#> 2   2.5  2.85
#> 3   3.5  3.15
```
