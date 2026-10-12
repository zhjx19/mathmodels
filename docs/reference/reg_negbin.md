# Negative Binomial Regression

Fits a negative binomial regression model (for over-dispersed counts)
via MASS::glm.nb().

## Usage

``` r
reg_negbin(formula, data, step = FALSE, ...)
```

## Arguments

- formula:

  A formula, e.g. `count ~ x1 + x2`.

- data:

  A data frame.

- step:

  Logical. Always `FALSE` for negative binomial.

- ...:

  Additional arguments passed to
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Value

A named list analogous to
[`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
plus:

- theta:

  Estimated dispersion parameter for the negative binomial distribution.

- residuals:

  Tibble of observed, fitted, Pearson residual.

## Details

**Note:** The `glm.nb` family does not support stepwise selection. Set
`step = FALSE`.

Unlike
[`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
the negative binomial model handles overdispersion through the estimated
dispersion parameter `theta` (fitted by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)): a small
`theta` indicates strong overdispersion, while `theta` growing large
approaches the Poisson model.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
mu = exp(1 + 0.3*x1)
y = rnbinom(100, size = 1, mu = mu)
res = reg_negbin(y ~ x1, data = data.frame(y, x1))
res$coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.907     0.117      7.74 1.01e-14   0.674      1.14 
#> 2 x1             0.322     0.116      2.78 5.52e- 3   0.0918     0.553
```
