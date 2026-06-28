# Back-Transform Forecasts to the Original Scale

Reverses the operations applied by
[`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
to restore point forecasts and confidence intervals to the original
scale. Differencing is un-done by cumulative summation using the tail of
the original series as the starting level.

## Usage

``` r
ts_back_transform(forecast_tbl, params, x_original)
```

## Arguments

- forecast_tbl:

  Tibble from
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md)
  (columns: `step`, `forecast`, optionally `lo_*` and `hi_*`).

- params:

  The `$params` element from a
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
  result.

- x_original:

  The original (pre-transform) `ts` or numeric vector, needed to supply
  the last observed level(s) for un-differencing.

## Value

The `forecast_tbl` with `forecast` and all interval columns converted to
the original scale.

## Examples

``` r
data(AirPassengers)
tr  = ts_transform(AirPassengers, method = "log", diff = 1,
                    seasonal_diff = 1, verbose = FALSE)
fit = ts_sarima(tr$transformed)
fc  = ts_forecast(fit, h = 12)
fc_orig = ts_back_transform(fc, tr$params, AirPassengers)
```
