# Two-Stage SARIMA-GARCH Joint Model

Stage 1: fits a SARIMA model on the level series via
[`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md).
Stage 2: fits a GARCH model on the SARIMA residuals via
[`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md).
Returns both sub-models plus a combined fitted tibble.

## Usage

``` r
ts_sarima_garch(
  x,
  sarima_order = NULL,
  sarima_seasonal = NULL,
  garch_order = c(1, 1),
  garch_dist = "std",
  frequency = NULL,
  ...
)
```

## Arguments

- x:

  A numeric vector or `ts` object.

- sarima_order:

  Integer vector `c(p,d,q)` or `NULL` (auto).

- sarima_seasonal:

  List `list(order=c(P,D,Q), period=m)` or `NULL`.

- garch_order:

  Integer vector `c(p,q)` (default `c(1,1)`).

- garch_dist:

  Character. Innovation distribution for GARCH (default `"std"`).

- frequency:

  Integer. Required for plain numeric `x`.

- ...:

  Additional arguments forwarded to
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md).

## Value

A named list:

- mean_model:

  Result list from
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md).

- variance_model:

  Result list from
  [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md).

- fitted:

  Combined tibble:
  `index, observed, mean_fitted, sigma, variance, std_residual`.

- model_info:

  Combined one-row summary tibble.

## Examples

``` r
data(AirPassengers)
res = ts_sarima_garch(log(AirPassengers))
res$model_info
#> # A tibble: 1 × 5
#>   mean_model            variance_model garch_dist sarima_aic garch_aic
#>   <chr>                 <chr>          <chr>           <dbl>     <dbl>
#> 1 ARIMA0,1,1(2,1,1)[12] GARCH(1,1)     std             -480.     -3.91
```
