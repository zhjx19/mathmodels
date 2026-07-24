# Complete Missing Time Points

Adds missing time points on a regular sequence and leaves their values
as `NA`.

## Usage

``` r
complete_ts_df(x, by)
```

## Arguments

- x:

  A `ts_df`.

- by:

  Step size passed to [`seq()`](https://rdrr.io/r/base/seq.html).

## Value

A completed `ts_df`.

## Examples

``` r
x = ts_df(data.frame(time = c(1, 2, 4), value = c(10, 12, 18)))
complete_ts_df(x, by = 1)
#>   time value
#> 1    1    10
#> 2    2    12
#> 3    3    NA
#> 4    4    18
```
