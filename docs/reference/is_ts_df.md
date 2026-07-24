# Check Whether an Object Is a `ts_df`

Check Whether an Object Is a `ts_df`

## Usage

``` r
is_ts_df(x)
```

## Arguments

- x:

  Object to check.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
x = ts_df(data.frame(time = 1:3, value = c(10, 12, 15)))
is_ts_df(x)
#> [1] TRUE
```
