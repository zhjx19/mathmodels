# Wilcoxon / Mann-Whitney Tests

Nonparametric counterparts of the t-test: Wilcoxon signed-rank
(one-sample and paired) and Mann-Whitney U (two-sample). Input shapes,
`.by` slicing and output layout match
[`stat_t_test()`](https://zhjx19.github.io/mathmodels/reference/stat_t_test.md).

## Usage

``` r
stat_wilcox_test(
  data,
  .cols,
  group,
  .by = NULL,
  mu = 0,
  paired = FALSE,
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
  Numeric test variable(s). One column for one-/two-sample tests;
  exactly two columns when `paired = TRUE`.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional grouping column with exactly two levels (the two groups to
  compare). Required for two-sample tests, forbidden for paired and
  one-sample tests.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; the test is repeated separately within each
  slice.

- mu:

  Null value of the mean (or mean difference).

- paired:

  Logical; paired test. Requires exactly two `.cols` and no `group`
  (differences are taken row-wise).

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- conf_level:

  Confidence level.

- p_adjust:

  P-value adjustment (see
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html)) applied
  across the returned rows.

## Value

A tibble with class `stat_infer` in the same layout as
[`stat_t_test()`](https://zhjx19.github.io/mathmodels/reference/stat_t_test.md),
with `statistic`, `p.value`, `effect` (rank-biserial correlation) and no
`df` column (NA). For one-sample tests the effect is the rank-biserial
correlation of the signed ranks versus `mu` (positive proportion minus
negative proportion); it is `NA` when all deviations from `mu` are zero.

## Examples

``` r
stat_wilcox_test(sleep, .cols = extra, group = group)
#> Wilcoxon rank-sum / signed-rank test (two.sided) 
#> # A tibble: 1 × 16
#>   variable n_used n_dropped estimate estimate1 estimate2 statistic    df p.value
#> * <chr>     <int>     <int>    <dbl>     <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 extra        20         0       NA        NA        NA      25.5    NA  0.0658
#> # ℹ 7 more variables: conf_low <dbl>, conf_high <dbl>, effect <dbl>,
#> #   method <chr>, alternative <chr>, p.adjusted <dbl>, sig <chr>
stat_wilcox_test(mtcars, .cols = mpg, mu = 20)
#> Wilcoxon rank-sum / signed-rank test (two.sided) 
#> # A tibble: 1 × 16
#>   variable n_used n_dropped estimate estimate1 estimate2 statistic    df p.value
#> * <chr>     <int>     <int>    <dbl>     <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 mpg          32         0       NA        NA        NA       249    NA   0.786
#> # ℹ 7 more variables: conf_low <dbl>, conf_high <dbl>, effect <dbl>,
#> #   method <chr>, alternative <chr>, p.adjusted <dbl>, sig <chr>
```
