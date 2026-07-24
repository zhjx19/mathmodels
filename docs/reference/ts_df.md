# Lightweight Time-Series Data Frame

Creates a simple single-series time-series data frame with fixed columns
`time` and `value`, plus `frequency` and `start` attributes.

## Usage

``` r
ts_df(data, time = "time", value = "value", frequency = 1, start = NULL)
```

## Arguments

- data:

  A data frame.

- time:

  Name of the time/index column.

- value:

  Name of the numeric value column.

- frequency:

  Positive numeric scalar. Defaults to 1.

- start:

  Optional start metadata. If `NULL`, uses the first sorted time.

## Value

A `ts_df` object, which is also a data frame.

## Examples

``` r
df = data.frame(time = c(2, 1, 3), value = c(12, 10, 15))
ts_df(df)
#>   time value
#> 1    1    10
#> 2    2    12
#> 3    3    15
```
