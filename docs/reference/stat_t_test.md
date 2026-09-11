# t-test (one-sample, two-sample, Welch, paired)

Performs t-tests on numeric columns of a tidy data frame. Each test
returns one row; multiple columns and `.by` slices are handled by row
repetition, with an optional p-value adjustment across rows.

## Usage

``` r
stat_t_test(
  data,
  .cols,
  group,
  .by = NULL,
  mu = 0,
  paired = FALSE,
  var.equal = FALSE,
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

- var.equal:

  Logical; pooled-variance t-test (`FALSE` gives Welch).

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- conf_level:

  Confidence level.

- p_adjust:

  P-value adjustment (see
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html)) applied
  across the returned rows.

## Value

A tibble with class `stat_infer`: `variable`, `.by` identifiers,
`n_used`, `n_dropped`, `estimate`, `estimate1`, `estimate2`,
`statistic`, `df`, `p.value`, `p.adjusted`, `sig`, `conf_low`,
`conf_high`, `effect` (Cohen's d), `method`, `alternative`. For
two-sample tests all difference quantities (`estimate`, `statistic`, CI,
`effect`) follow the group1 - group2 direction (first level minus second
level, rstatix convention).

## Examples

``` r
stat_t_test(sleep, .cols = extra, group = group)
#> t-test (two.sided) 
#> # A tibble: 1 × 16
#>   variable n_used n_dropped estimate estimate1 estimate2 statistic    df p.value
#> * <chr>     <int>     <int>    <dbl>     <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 extra        20         0    -1.58      0.75      2.33     -1.86  17.8  0.0794
#> # ℹ 7 more variables: conf_low <dbl>, conf_high <dbl>, effect <dbl>,
#> #   method <chr>, alternative <chr>, p.adjusted <dbl>, sig <chr>
stat_t_test(mtcars, .cols = mpg, mu = 20)
#> t-test (two.sided) 
#> # A tibble: 1 × 16
#>   variable n_used n_dropped estimate estimate1 estimate2 statistic    df p.value
#> * <chr>     <int>     <int>    <dbl>     <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 mpg          32         0     20.1        NA        NA    0.0851    31   0.933
#> # ℹ 7 more variables: conf_low <dbl>, conf_high <dbl>, effect <dbl>,
#> #   method <chr>, alternative <chr>, p.adjusted <dbl>, sig <chr>
stat_t_test(mtcars, .cols = mpg, group = am, var.equal = TRUE)
#> t-test (two.sided) 
#> # A tibble: 1 × 16
#>   variable n_used n_dropped estimate estimate1 estimate2 statistic    df p.value
#> * <chr>     <int>     <int>    <dbl>     <dbl>     <dbl>     <dbl> <dbl>   <dbl>
#> 1 mpg          32         0    -7.24      17.1      24.4     -4.11    30 2.85e-4
#> # ℹ 7 more variables: conf_low <dbl>, conf_high <dbl>, effect <dbl>,
#> #   method <chr>, alternative <chr>, p.adjusted <dbl>, sig <chr>
```
