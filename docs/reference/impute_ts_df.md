# Complete and Interpolate a `ts_df`

Completes missing time points and fills internal missing values by
linear or spline interpolation. Endpoint values must already be
observed.

## Usage

``` r
impute_ts_df(x, by, method = c("linear", "spline"))
```

## Arguments

- x:

  A `ts_df`.

- by:

  Step size passed to [`seq()`](https://rdrr.io/r/base/seq.html).

- method:

  Interpolation method: `"linear"` or `"spline"`.

## Value

An interpolated `ts_df`.

## Examples

``` r
x = ts_df(data.frame(time = c(1, 3, 5), value = c(10, 20, 30)))
impute_ts_df(x, by = 1, method = "linear")
#>   time value
#> 1    1    10
#> 2    2    15
#> 3    3    20
#> 4    4    25
#> 5    5    30
```
