# Linear Interpolation

Performs linear (piecewise) interpolation at specified points.

## Usage

``` r
interp_linear(x, y, xout)
```

## Arguments

- x:

  Numeric vector of observed x-values (no duplicates).

- y:

  Numeric vector of observed y-values, same length as `x`.

- xout:

  Numeric vector of x-values where interpolation is desired.

## Value

A tibble with columns `x` and `y`, containing the interpolated values at
`xout`. The attribute `method` is set to `"linear"`.

## Examples

``` r
x = c(1, 2, 3, 4, 5)
y = c(2, 4, 1, 5, 3)
interp_linear(x, y, xout = c(1.5, 2.5, 3.5))
#> # A tibble: 3 × 2
#>       x     y
#>   <dbl> <dbl>
#> 1   1.5   3  
#> 2   2.5   2.5
#> 3   3.5   3  
```
