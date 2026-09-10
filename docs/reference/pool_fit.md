# Pool lm/glm Fits Across Imputations

Extracts coefficients and standard errors from a list of `lm`/`glm` fits
(one per imputation, usually built with
[`impute_multiple()`](https://zhjx19.github.io/mathmodels/reference/impute_multiple.md))
and pools each coefficient via
[`rubin_pool()`](https://zhjx19.github.io/mathmodels/reference/rubin_pool.md),
automatically passing the complete-data residual degrees of freedom so
the Barnard-Rubin small-sample adjustment applies.

## Usage

``` r
pool_fit(fits, alpha = 0.05)
```

## Arguments

- fits:

  A list of `lm`/`glm` fits, same model on the imputed data.

- alpha:

  Significance level. Default 0.05.

## Value

A tibble with one row per coefficient:
`term, estimate, std_error, lcl, ucl, df, t, p_value, fmi, r, rel_eff`.

## Examples

``` r
set.seed(1)
d = data.frame(x = rnorm(40), z = rnorm(40))
d$y = 0.5 * d$x + rnorm(40)
d$x[c(4, 9)] = NA
imps = impute_multiple(d, .cols = x, m = 5, maxit = 3, seed = 2)
fits = lapply(imps, \(dat) lm(y ~ x, dat))
pool_fit(fits)
#> Pooled estimates (Rubin's rules)
#> # A tibble: 2 × 7
#>   term        estimate std_error     lcl   ucl p_value    fmi
#>   <chr>          <dbl>     <dbl>   <dbl> <dbl>   <dbl>  <dbl>
#> 1 (Intercept)    0.151     0.130 -0.113  0.415  0.255  0.0565
#> 2 x              0.327     0.153  0.0155 0.638  0.0402 0.0557
```
