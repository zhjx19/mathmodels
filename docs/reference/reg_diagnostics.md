# Model Diagnostics

Runs diagnostic tests on a fitted regression model.

## Usage

``` r
reg_diagnostics(model_result)
```

## Arguments

- model_result:

  A list returned by any `reg_*()` function.

## Value

A named list with diagnostic tibbles. For linear models: VIF,
Shapiro-Wilk normality, Breusch-Pagan, Durbin-Watson, residual stats.
For logistic/poisson/negative binomial: Hosmer-Lemeshow / dispersion /
deviance stats.
