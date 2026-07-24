# Piecewise Cubic Hermite Interpolation

Shape-preserving (monotone) piecewise cubic Hermite interpolation using
the Fritsch-Carlson method via `splinefun(method = "monoH.FC")`.

## Usage

``` r
interp_hermite(x, y, xout)
```

## Arguments

- x:

  Numeric vector of observed x-values (no duplicates, length \>= 2).

- y:

  Numeric vector of observed y-values.

- xout:

  Numeric vector of x-values where interpolation is desired.

## Value

A tibble with columns `x` and `y`. The attribute `method` is set to
`"hermite"`.

## Examples

``` r
x = c(1, 2, 3, 4, 5)
y = c(1, 2, 4, 7, 11)
interp_hermite(x, y, xout = c(1.5, 2.5, 3.5))
#> # A tibble: 3 × 2
#>       x     y
#>   <dbl> <dbl>
#> 1   1.5  1.44
#> 2   2.5  2.88
#> 3   3.5  5.38
```
