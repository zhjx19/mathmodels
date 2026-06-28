# Residual Diagnostic Plots

Four-panel diagnostic: (1) residuals over time, (2) ACF of residuals,
(3) Q-Q plot, (4) histogram vs normal density.

## Usage

``` r
plot_ts_residuals(model_result, title = "Residual Diagnostics", max_lag = 30L)
```

## Arguments

- model_result:

  Result from any `ts_*()` function that contains a `$fitted` tibble
  with a `residual` column.

- title:

  Character. Plot title.

- max_lag:

  Integer. Maximum lag for ACF (default `30`).

## Value

A `patchwork` composite `ggplot`.

## Examples

``` r
data(AirPassengers)
fit = ts_sarima(log(AirPassengers))
plot_ts_residuals(fit)

```
