# ETS (Error, Trend, Seasonality) Exponential Smoothing

Wraps
[`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
and returns a tidy result list.

## Usage

``` r
ts_ets(x, model = "ZZZ", frequency = NULL, ...)
```

## Arguments

- x:

  A numeric vector or `ts` object.

- model:

  Character ETS model string, e.g. `"ZZZ"` (auto-select, default),
  `"AAN"`, `"AAA"`.

- frequency:

  Integer. Required when `x` is a plain numeric vector.

- ...:

  Additional arguments forwarded to
  [`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html).

## Value

A named list:

- model_info:

  One-row tibble: model type, AIC, AICc, BIC, log-likelihood.

- parameters:

  Tidy tibble of smoothing parameters and initial states.

- fitted:

  Tibble with `index, observed, fitted, residual`.

- model:

  Raw `ets` object.

## Examples

``` r
data(AirPassengers)
res = ts_ets(AirPassengers)
res$model_info
#> # A tibble: 1 × 5
#>   model_type  log_lik   aic  aicc   bic
#>   <chr>         <dbl> <dbl> <dbl> <dbl>
#> 1 ETS(M,Ad,M)   -680. 1395. 1401. 1449.
```
