# Drop Missing Values from a `ts_df`

Drop Missing Values from a `ts_df`

## Usage

``` r
drop_na_ts_df(x)
```

## Arguments

- x:

  A `ts_df`.

## Value

A `ts_df` with rows whose `value` is `NA` removed.

## Examples

``` r
x = ts_df(data.frame(time = 1:4, value = c(10, NA, 12, 13)))
drop_na_ts_df(x)
#>   time value
#> 1    1    10
#> 2    3    12
#> 3    4    13
```
