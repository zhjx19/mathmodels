# Logistic Regression

Fits a binary logistic regression model via `glm(family = binomial())`,
with optional stepwise selection.

## Usage

``` r
reg_logistic(
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

- odds_ratio:

  Tibble with odds ratios and 95% CIs.

- residuals:

  Tibble of observed (0/1), fitted probability, Pearson residual.

## Details

Same stepwise logic as
[`reg_lm`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md), but
uses GLM with the binomial family. A Hosmer-Lemeshow goodness-of-fit
test is included in diagnostics.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
logit = 0.5 + x1
prob = 1 / (1 + exp(-logit))
y = rbinom(100, 1, prob)
res = reg_logistic(y ~ x1, data = data.frame(y, x1))
res$coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic   p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.285     0.226      1.26 0.207       -0.163     0.734
#> 2 x1             1.07      0.270      3.98 0.0000682    0.539     1.61 
```
