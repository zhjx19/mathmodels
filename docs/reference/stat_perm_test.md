# Permutation Test

Generic permutation-test framework. With `statistic = NULL` the default
statistic is chosen by `type`: group mean difference (`two.sample`),
mean of differences (`paired`), mean around `mu` (`one.sample`), or
Pearson correlation (`independence`). A custom `statistic` of the form
`function(data, indices, ...)` is supported for `two.sample` (indices =
randomly re-ordered row indices, free permutation) and `independence`
(indices = the re-ordered rows used to shuffle the second column).

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
  for one-/two-sample, two columns for paired/independence.

- group:

  Two-level factor column (required for `two.sample` when
  `statistic = NULL`).

- statistic:

  Optional custom function `function(data, indices, ...)` returning a
  scalar; observed value uses `indices = seq_len(n)`.

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
