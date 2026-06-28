# SARIMA Model Fitting

Auto-selects or fits a user-specified SARIMA model via
[`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
or
[`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html),
returning a comprehensive tidy result.

## Usage

``` r
ts_sarima(
  x,
  order = NULL,
  seasonal = NULL,
  frequency = NULL,
  stepwise = TRUE,
  approximation = TRUE,
  ...
)
```

## Arguments

- x:

  A numeric vector or `ts` object.

- order:

  Integer vector of length 3: `c(p, d, q)`. `NULL` (default) triggers
  automatic selection.

- seasonal:

  A list `list(order = c(P,D,Q), period = m)` or `NULL` for automatic
  selection.

- frequency:

  Integer. Required when `x` is a plain vector.

- stepwise:

  Logical. Use stepwise search in auto.arima (default `TRUE`).

- approximation:

  Logical. Use approximation in auto.arima (default `TRUE`).

- ...:

  Additional arguments forwarded to
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
  or
  [`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html).

## Value

A named list:

- model_info:

  One-row tibble: ARIMA order string, AIC, AICc, BIC, log-likelihood,
  sigma².

- coefficients:

  Tidy tibble with estimate and std error.

- fitted:

  Tibble: `index, observed, fitted, residual`.

- diagnostics:

  Ljung-Box test result tibble.

- model:

  Raw `Arima` object.

## Examples

``` r
data(AirPassengers)
res = ts_sarima(log(AirPassengers))
res$model_info
#> # A tibble: 1 × 6
#>   model_type            log_lik   aic  aicc   bic  sigma2
#>   <chr>                   <dbl> <dbl> <dbl> <dbl>   <dbl>
#> 1 ARIMA0,1,1(2,1,1)[12]    245. -480. -479. -466. 0.00138
```
