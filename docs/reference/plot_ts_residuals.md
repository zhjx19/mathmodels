# Residual Diagnostic Plot

Residual Diagnostic Plot

## Usage

``` r
plot_ts_residuals(model_result, max_lag = 30L, title = "Residual Diagnostics")
```

## Arguments

- model_result:

  Result from
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md)
  or
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md).

- max_lag:

  Maximum lag.

- title:

  Plot title.

## Value

A patchwork plot.

## Examples

``` r
fit = ts_sarima(as_ts_df(log(AirPassengers)), stepwise = TRUE, approximation = TRUE)
plot_ts_residuals(fit)
```
