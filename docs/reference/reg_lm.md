# Multivariable Linear Regression

Fits an OLS model via [`lm()`](https://rdrr.io/r/stats/lm.html),
optionally performing stepwise selection using Akaike's information
criterion.

## Usage

``` r
reg_lm(
  formula,
  data,
  step = FALSE,
  direction = c("both", "forward", "backward"),
  ...
)
```

## Arguments

- formula:

  A formula, e.g. `y ~ x1 + x2`.

- data:

  A data frame containing the variables in `formula`.

- step:

  Logical. Whether to perform stepwise selection (`TRUE`/`FALSE`).
  Default is `FALSE`.

- direction:

  Character: `"both"` (default), `"forward"`, or `"backward"`.

- ...:

  Additional arguments passed to
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`step()`](https://rdrr.io/r/stats/step.html).

## Value

A named list:

- model:

  Fitted `lm` object.

- coefficient:

  Tibble of coefficients, standard errors, t-values, p-values, 95% CIs.

- model_info:

  Tibble of fit statistics (`r.squared`, `adj.r.squared`, AIC, BIC).

- residuals:

  Tibble of observed, fitted, residuals.

- diagnostics:

  Tibble of key fit metrics (`Breusch-Pagan`, Durbin-Watson).

- formula:

  Original formula.

- input:

  The original data frame.

## Details

When `step = TRUE`, the function first fits a null model (response
only), then uses [`step()`](https://rdrr.io/r/stats/step.html) to search
for the best combination of terms by AIC.

## Examples

``` r
set.seed(42)
x1 = rnorm(100); x2 = rnorm(100)
y = 1 + 2*x1 - 3*x2 + rnorm(100, sd = 2)
res = reg_lm(y ~ x1 + x2, data = data.frame(y, x1, x2))
res$coefficient
#> # A tibble: 3 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)     1.00     0.204      4.92 3.48e- 6    0.599      1.41
#> 2 x1              1.71     0.196      8.75 6.72e-14    1.32       2.10
#> 3 x2             -2.83     0.225    -12.5  4.91e-22   -3.28      -2.38
```
