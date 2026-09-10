# Correlation Analysis

Pairwise correlations among numeric columns with test statistics,
p-values and (Pearson) confidence intervals. Each pair uses its own
pairwise-complete rows (rows where both variables are observed, the
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)
convention); `n_used` is the per-pair complete count.

## Usage

``` r
stat_cor(
  data,
  .cols,
  method = c("pearson", "spearman", "kendall"),
  .by = NULL,
  alternative = "two.sided",
  conf_level = 0.95,
  p_adjust = "none"
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns (at least 2).

- method:

  `"pearson"`, `"spearman"`, or `"kendall"`. Confidence intervals are
  only available for Pearson (`NA` otherwise).

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- conf_level:

  Confidence level (Pearson only).

- p_adjust:

  P-value adjustment across returned rows.

## Value

A tibble with class `stat_infer`: `.by` identifiers, `var1`, `var2`,
`n_used`, `estimate`, `statistic`, `df`, `p.value`, `p.adjusted`, `sig`,
`conf_low`, `conf_high`, `method`. Each row is computed on the
pairwise-complete rows of its two variables (`n_used` = that pair's
complete-row count).

## Examples

``` r
stat_cor(mtcars, .cols = mpg:wt)
#> Correlation test (pearson, two.sided) 
#> # A tibble: 15 × 12
#>    var1  var2  n_used estimate statistic    df  p.value conf_low conf_high
#>  * <chr> <chr>  <int>    <dbl>     <dbl> <int>    <dbl>    <dbl>     <dbl>
#>  1 mpg   cyl       32   -0.852     -8.92    30 6.11e-10   -0.926    -0.716
#>  2 mpg   disp      32   -0.848     -8.75    30 9.38e-10   -0.923    -0.708
#>  3 mpg   hp        32   -0.776     -6.74    30 1.79e- 7   -0.885    -0.586
#>  4 mpg   drat      32    0.681      5.10    30 1.78e- 5    0.436     0.832
#>  5 mpg   wt        32   -0.868     -9.56    30 1.29e-10   -0.934    -0.744
#>  6 cyl   disp      32    0.902     11.4     30 1.80e-12    0.807     0.951
#>  7 cyl   hp        32    0.832      8.23    30 3.48e- 9    0.682     0.915
#>  8 cyl   drat      32   -0.700     -5.37    30 8.24e- 6   -0.843    -0.465
#>  9 cyl   wt        32    0.782      6.88    30 1.22e- 7    0.597     0.889
#> 10 disp  hp        32    0.791      7.08    30 7.14e- 8    0.611     0.893
#> 11 disp  drat      32   -0.710     -5.53    30 5.28e- 6   -0.849    -0.481
#> 12 disp  wt        32    0.888     10.6     30 1.22e-11    0.781     0.944
#> 13 hp    drat      32   -0.449     -2.75    30 9.99e- 3   -0.690    -0.119
#> 14 hp    wt        32    0.659      4.80    30 4.15e- 5    0.403     0.819
#> 15 drat  wt        32   -0.712     -5.56    30 4.78e- 6   -0.850    -0.484
#> # ℹ 3 more variables: method <chr>, p.adjusted <dbl>, sig <chr>
stat_cor(mtcars, .cols = mpg:wt, method = "spearman")
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Correlation test (spearman, two.sided) 
#> # A tibble: 15 × 12
#>    var1  var2  n_used estimate statistic    df  p.value conf_low conf_high
#>  * <chr> <chr>  <int>    <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>
#>  1 mpg   cyl       32   -0.911    10425.    NA 4.69e-13       NA        NA
#>  2 mpg   disp      32   -0.909    10415.    NA 6.37e-13       NA        NA
#>  3 mpg   hp        32   -0.895    10337.    NA 5.09e-12       NA        NA
#>  4 mpg   drat      32    0.651     1902.    NA 5.38e- 5       NA        NA
#>  5 mpg   wt        32   -0.886    10292.    NA 1.49e-11       NA        NA
#>  6 cyl   disp      32    0.928      395.    NA 2.28e-14       NA        NA
#>  7 cyl   hp        32    0.902      536.    NA 1.87e-12       NA        NA
#>  8 cyl   drat      32   -0.679     9160.    NA 1.94e- 5       NA        NA
#>  9 cyl   wt        32    0.858      776.    NA 3.57e-10       NA        NA
#> 10 disp  hp        32    0.851      813.    NA 6.79e-10       NA        NA
#> 11 disp  drat      32   -0.684     9186.    NA 1.61e- 5       NA        NA
#> 12 disp  wt        32    0.898      558.    NA 3.35e-12       NA        NA
#> 13 hp    drat      32   -0.520     8294.    NA 2.28e- 3       NA        NA
#> 14 hp    wt        32    0.775     1229.    NA 1.95e- 7       NA        NA
#> 15 drat  wt        32   -0.750     9550.    NA 7.59e- 7       NA        NA
#> # ℹ 3 more variables: method <chr>, p.adjusted <dbl>, sig <chr>
```
