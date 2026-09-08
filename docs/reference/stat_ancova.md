# Analysis of Covariance (ANCOVA)

Fits `resp ~ group + covariates` (group as factor, covariates numeric)
and returns the ANOVA table with partial eta-squared. Optionally checks
the parallel-slopes assumption (group-by-covariate interactions).

## Usage

``` r
stat_ancova(data, resp, group, covariates, .by = NULL, test_parallel = FALSE)
```

## Arguments

- data:

  A data frame.

- resp:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric response column(s); each column is analysed separately.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Single factor column.

- covariates:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric covariate column(s).

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

- test_parallel:

  Logical; if `TRUE`, the group-by-covariate interactions are tested and
  a warning is issued when any is significant at 0.05 (parallel-slopes
  assumption in doubt).

## Value

A tibble with class `stat_infer`: `resp`, `.by` identifiers, `term`,
`df`, `sumsq`, `meansq`, `statistic`, `p.value`, `sig`, `effect`
(partial eta-squared), `method`.

## Examples

``` r
d = data.frame(y = c(1:10, 2:11), g = rep(c("a", "b"), each = 10),
               x = 1:20)
stat_ancova(d, resp = y, group = g, covariates = x)
#> Warning: ANOVA F-tests on an essentially perfect fit are unreliable
#> Analysis of covariance (ANCOVA) 
#> # A tibble: 3 × 9
#>   resp  term         df    sumsq   meansq statistic    p.value effect sig  
#> * <chr> <chr>     <int>    <dbl>    <dbl>     <dbl>      <dbl>  <dbl> <chr>
#> 1 y     g             1 5.00e+ 0 5.00e+ 0   1.66e30  7.31e-248      1 "***"
#> 2 y     x             1 1.65e+ 2 1.65e+ 2   5.48e31  9.05e-261      1 "***"
#> 3 y     Residuals    17 5.11e-29 3.01e-30  NA       NA             NA ""   
```
