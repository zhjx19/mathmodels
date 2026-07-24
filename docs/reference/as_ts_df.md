# Convert Common Inputs to `ts_df`

Converts a `ts` object, numeric vector, explicit data-frame columns, or
an existing `ts_df` to the package's lightweight time-series structure.

## Usage

``` r
as_ts_df(x, time = NULL, value = NULL, frequency = NULL, start = NULL)
```

## Arguments

- x:

  A `ts`, numeric vector, data frame, or `ts_df`.

- time, value:

  Required column names when `x` is a data frame.

- frequency:

  Optional frequency metadata.

- start:

  Optional start metadata.

## Value

A `ts_df` object.

## Examples

``` r
as_ts_df(ts(1:4, frequency = 4))
#>   time value
#> 1    1     1
#> 2    2     2
#> 3    3     3
#> 4    4     4
as_ts_df(c(3, 5, 8))
#>   time value
#> 1    1     3
#> 2    2     5
#> 3    3     8
as_ts_df(data.frame(month = 1:3, demand = c(10, 12, 15)), time = "month", value = "demand")
#>   time value
#> 1    1    10
#> 2    2    12
#> 3    3    15
```
