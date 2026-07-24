# Back-Transform Forecasts

Converts forecasts produced on a transformed scale back to the original
scale.

## Usage

``` r
ts_back_transform(forecast_tbl, params)
```

## Arguments

- forecast_tbl:

  Forecast tibble from
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md).

- params:

  The `params` element from
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md).

## Value

The forecast tibble with numeric forecast columns restored.

## Examples

``` r
x = as_ts_df(AirPassengers)
tr = ts_transform(x, method = "log")
fc = tibble::tibble(step = 1:2, forecast = log(c(500, 600)))
ts_back_transform(fc, tr$params)
#> # A tibble: 2 × 2
#>    step forecast
#>   <int>    <dbl>
#> 1     1      500
#> 2     2      600
```
