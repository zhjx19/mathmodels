# Poisson Regression

Fits a count-response regression model via `glm(family = poisson())`,
with optional stepwise selection.

## Usage

``` r
reg_poisson(
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

A named list analogous to
[`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
plus:

- dispersion:

  Dispersion ratio (should be near 1).

- residuals:

  Tibble with observed, fitted, Pearson residual.

## Details

Performs an overdispersion check: dispersion \> 1.5 suggests NB may be
preferable.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
mu = exp(1 + 0.3*x1)
y = rpois(100, mu)
res = reg_poisson(y ~ x1, data = data.frame(y, x1))
res$coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.848    0.0687     12.3  5.32e-35    0.712     0.985
#> 2 x1             0.423    0.0664      6.37 1.92e-10    0.291     0.554
```
