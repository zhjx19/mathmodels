# Generate Forecasts

Generate Forecasts

## Usage

``` r
ts_forecast(model_result, h = 12, level = c(80, 95), newxreg = NULL)
```

## Arguments

- model_result:

  Result from
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
  or
  [`ts_arimax()`](https://zhjx19.github.io/mathmodels/reference/ts_arimax.md).

- h:

  Forecast horizon.

- level:

  Confidence levels.

- newxreg:

  Future exogenous regressors for ARIMAX models. A matrix or data.frame
  with `nrow(newxreg) == h` and the same number of columns as the
  original `xreg`.

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

x = as_ts_df(log(AirPassengers))
xreg = data.frame(trend = seq_len(nrow(x)))
fit_arimax = ts_arimax(x, xreg = xreg, order = c(1, 1, 1))
newxreg = data.frame(trend = nrow(x) + 1:3)
ts_forecast(fit_arimax, h = 3, newxreg = newxreg)
#> # A tibble: 3 × 6
#>    step forecast lo_80 hi_80 lo_95 hi_95
#>   <int>    <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1     1     6.15  6.02  6.28  5.95  6.35
#> 2     2     6.12  5.91  6.33  5.79  6.44
#> 3     3     6.15  5.89  6.41  5.76  6.54
```
