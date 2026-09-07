# Back-Transform Forecasts

Converts forecasts produced on a transformed scale back to the original
scale.

## Usage

``` r
ts_back_transform(forecast_tbl, params, bias_adj = TRUE)
```

## Arguments

- forecast_tbl:

  Forecast tibble from
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md).

- params:

  The `params` element from
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md).

- bias_adj:

  Logical; if `TRUE` (default), the point-forecast column is
  back-transformed with a bias adjustment so that forecasts on the
  log/Box-Cox scale represent means rather than medians. The forecast
  variance required by the adjustment is derived per horizon from the
  widest stored prediction interval, so it applies only to genuine
  forecast rows: when the tibble has no `lo_*`/`hi_*` interval columns
  the plain transform is used (historical/fitted values round-trip
  unchanged). Interval columns are quantiles and are always transformed
  with the plain transform.

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
