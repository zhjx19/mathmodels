# Generate Forecasts

Generate Forecasts

## Usage

``` r
ts_forecast(model_result, h = 12, level = c(80, 95))
```

## Arguments

- model_result:

  Result from
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md)
  or
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md).

- h:

  Forecast horizon.

- level:

  Confidence levels.

## Value

A tibble with point forecasts and intervals.

## Examples

``` r
fit = ts_ets(as_ts_df(log(AirPassengers)))
ts_forecast(fit, h = 3)
#> # A tibble: 3 × 6
#>    step forecast lo_80 hi_80 lo_95 hi_95
#>   <int>    <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1     1     6.11  6.05  6.16  6.02  6.19
#> 2     2     6.10  6.03  6.16  6.00  6.19
#> 3     3     6.25  6.18  6.32  6.14  6.36
```
