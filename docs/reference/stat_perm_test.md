# Permutation Test

Generic permutation-test framework. With `statistic = NULL` the default
statistic is chosen by `type`: group mean difference (`two.sample`),
mean of differences (`paired`), mean around `mu` (`one.sample`), or
Pearson correlation (`independence`). A custom `statistic` of the form
`function(data, indices, ...)` is supported for `two.sample` and
`independence`; see Details for the permutation semantics.

## Usage

``` r
stat_perm_test(
  data,
  .cols = NULL,
  group = NULL,
  statistic = NULL,
  ...,
  type = c("two.sample", "paired", "one.sample", "independence"),
  alternative = "two.sided",
  permutations = 2000,
  mu = 0,
  seed = NULL,
  .by = NULL
)
```

## Arguments

- data:

  A data frame.

- .cols:

  Numeric test column(s) (required when `statistic = NULL`): one column
  for one-/two-sample, two columns for paired/independence. With a
  custom statistic, `independence` requires exactly two `.cols` columns
  (the second one is shuffled).

- group:

  Two-level factor column. Required for `two.sample` (also with a custom
  statistic: its labels are permuted).

- statistic:

  Optional custom function `function(data, indices, ...)` returning a
  scalar; the observed value uses `indices = seq_len(n)`.

- ...:

  Additional arguments passed to `statistic`.

- type:

  `"two.sample"`, `"paired"`, `"one.sample"`, or `"independence"`.

- alternative:

  `"two.sided"`, `"greater"`, or `"less"`.

- permutations:

  Number of permutations.

- mu:

  Null value for `one.sample`.

- seed:

  Optional integer seed.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; the test is repeated within each slice.

## Value

A tibble with class `stat_infer`: `.by` identifiers, `observed`,
`p.value`, `null_mean`, `null_sd`, `permutations`, `type`,
`alternative`.

## Details

Custom-statistic semantics: for `type = "two.sample"` the null
distribution is built by permuting the GROUP LABELS — for every
permutation the `group` column of the data handed to `statistic` is
re-assigned at random (`indices = seq_len(n)`), so any group-contrast
statistic varies across permutations. For `type = "independence"` the
SECOND `.cols` column is shuffled (within each `.by` slice; slice
columns are never touched) before `statistic` sees the data. If the
permuted statistic is (near-)constant — for example a statistic that
does not depend on the permuted labels — a warning is issued instead of
silently reporting a meaningless p = 1.

## Examples

``` r
set.seed(1)
d = data.frame(y = c(rnorm(20), rnorm(20, 0.8)),
               g = rep(c("a", "b"), each = 20))
stat_perm_test(d, .cols = y, group = g, type = "two.sample")
#> Permutation test (two.sample, two.sided) 
#> # A tibble: 1 × 8
#>   observed p.value null_mean null_sd permutations type       alternative sig  
#> *    <dbl>   <dbl>     <dbl>   <dbl>        <dbl> <chr>      <chr>       <chr>
#> 1   -0.603  0.0365   0.00270   0.293         2000 two.sample two.sided   *    
```
