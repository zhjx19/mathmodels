# Polynomial Interpolation

Fits a polynomial of specified degree via
[`lm()`](https://rdrr.io/r/stats/lm.html) and evaluates it at `xout`.

## Usage

``` r
interp_poly(x, y, xout, degree = 3)
```

## Arguments

- x:

  Numeric vector of observed x-values (no duplicates).

- y:

  Numeric vector of observed y-values.

- xout:

  Numeric vector of x-values where interpolation is desired.

- degree:

  Integer, degree of the polynomial. Minimum 1.

## Value

A tibble with columns `x` and `y`. The attribute `method` is set to
`"poly"`.

## Examples

``` r
x = c(1, 2, 3, 4, 5)
y = c(2, 4, 1, 5, 3)
interp_poly(x, y, xout = c(1.5, 2.5, 3.5), degree = 3)
#> # A tibble: 3 × 2
#>       x     y
#>   <dbl> <dbl>
#> 1   1.5  2.39
#> 2   2.5  2.84
#> 3   3.5  3.41
```
