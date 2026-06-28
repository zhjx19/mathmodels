# GARCH Variance Modeling

Fits a GARCH(p,q) model (optionally with ARMA mean) using
[`rugarch::ugarchfit()`](https://rdrr.io/pkg/rugarch/man/ugarchfit-methods.html)
and returns tidy results.

## Usage

``` r
ts_garch(
  x,
  garch_order = c(1, 1),
  arma_order = c(0, 0),
  dist = "norm",
  frequency = NULL,
  ...
)
```

## Arguments

- x:

  A numeric vector or `ts` object (typically returns/residuals).

- garch_order:

  Integer vector `c(p, q)`: GARCH order (default `c(1, 1)`).

- arma_order:

  Integer vector `c(p, q)`: ARMA order for the mean equation (default
  `c(0, 0)` = constant mean).

- dist:

  Character. Innovation distribution: `"norm"` (default), `"std"`
  (Student-t), `"ged"`, `"snorm"`, `"sstd"`.

- frequency:

  Integer. Required when `x` is a plain vector.

- ...:

  Additional arguments forwarded to
  [`rugarch::ugarchspec()`](https://rdrr.io/pkg/rugarch/man/ugarchspec-methods.html).

## Value

A named list:

- model_info:

  One-row tibble: spec, distribution, log-likelihood, AIC, BIC.

- coefficients:

  Tidy coefficient tibble with estimates, std errors, t-stats, p-values.

- fitted:

  Tibble: `index, observed, sigma, variance, std_residual`.

- diagnostics:

  ARCH-LM and Ljung-Box tests on squared residuals.

- model:

  Raw `uGARCHfit` object.

## Examples

``` r
set.seed(42)
r = diff(log(AirPassengers))
res = ts_garch(r)
res$model_info
#> # A tibble: 1 × 5
#>   model_type distribution log_lik   aic   bic
#>   <chr>      <chr>          <dbl> <dbl> <dbl>
#> 1 GARCH(1,1) norm            118. -1.59 -1.51
```
