# Plot Forecasts

Plot Forecasts

## Usage

``` r
plot_ts_forecast(
  x,
  forecast_tbl,
  level = c(80, 95),
  by = NULL,
  title = "Forecast"
)
```

## Arguments

- x:

  Original `ts_df`.

- forecast_tbl:

  Result from
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md).

- level:

  Confidence levels to draw when present.

- by:

  Step size for future time values.

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
fit = ts_ets(as_ts_df(log(AirPassengers)))
fc = ts_forecast(fit, h = 3)
plot_ts_forecast(as_ts_df(log(AirPassengers)), fc, by = 1)
```
