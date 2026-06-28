# Generate Forecasts from a Fitted Time Series Model

A unified interface for forecasting from ETS, SARIMA, GARCH, or
SARIMA-GARCH model results produced by the `ts_*()` functions.

## Usage

``` r
ts_forecast(model_result, h = 12, level = c(80, 95), x_ts = NULL)
```

## Arguments

- model_result:

  A result list from
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
  [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md),
  or
  [`ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima_garch.md).

- h:

  Integer. Forecast horizon (number of steps ahead).

- level:

  Numeric vector. Confidence levels in percent (default `c(80, 95)`).

- x_ts:

  Original `ts` object (required for GARCH/SARIMA-GARCH to preserve time
  attributes). If `NULL`, integer steps are used.

## Value

A `tibble` with columns:

- step:

  Forecast step (1 to h).

- forecast:

  Point forecast.

- lo\_:

  Lower bound for each confidence level.

- hi\_:

  Upper bound for each confidence level.

- sigma:

  Predicted conditional std dev (GARCH models only).

## Examples

``` r
data(AirPassengers)
fit = ts_sarima(log(AirPassengers))
ts_forecast(fit, h = 24)
#> # A tibble: 24 × 6
#>     step forecast lo_80 hi_80 lo_95 hi_95
#>    <int>    <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1     1     6.11  6.06  6.16  6.04  6.18
#>  2     2     6.05  6.00  6.11  5.97  6.14
#>  3     3     6.17  6.11  6.24  6.08  6.27
#>  4     4     6.20  6.13  6.27  6.09  6.30
#>  5     5     6.23  6.16  6.31  6.12  6.35
#>  6     6     6.37  6.29  6.45  6.25  6.49
#>  7     7     6.51  6.43  6.59  6.38  6.64
#>  8     8     6.51  6.42  6.60  6.37  6.64
#>  9     9     6.33  6.23  6.42  6.18  6.47
#> 10    10     6.21  6.11  6.31  6.06  6.36
#> # ℹ 14 more rows
```
