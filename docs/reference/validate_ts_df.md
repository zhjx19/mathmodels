# Validate a `ts_df`

Validate a `ts_df`

## Usage

``` r
validate_ts_df(x, require_complete = TRUE, require_no_na = FALSE)
```

## Arguments

- x:

  A `ts_df` object.

- require_complete:

  If `TRUE`, first and last values must not be `NA`.

- require_no_na:

  If `TRUE`, no values may be `NA`.

## Value

A sorted `ts_df` object.

## Examples

``` r
x = ts_df(data.frame(time = 1:3, value = c(10, 12, 15)))
validate_ts_df(x)
#>   time value
#> 1    1    10
#> 2    2    12
#> 3    3    15
```
