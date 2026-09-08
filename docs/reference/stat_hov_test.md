# Variance Homogeneity Tests

Tests the equality of variances across the levels of `group`: Levene
(via [`car::leveneTest`](https://rdrr.io/pkg/car/man/leveneTest.html),
centering on median or mean), Bartlett (assumes normality), or the
two-level F test.

## Usage

``` r
stat_hov_test(
  data,
  .cols,
  group,
  .by = NULL,
  method = c("levene", "bartlett", "f"),
  center = c("median", "mean")
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric response column(s); each column is analysed separately.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Single factor column.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

- method:

  `"levene"`, `"bartlett"`, or `"f"`.

- center:

  For Levene: `"median"` (Brown-Forsythe, default) or `"mean"`.

## Value

A tibble with class `stat_infer`: `variable`, `.by` identifiers,
`n_used`, `statistic`, `df`, `df2`, `p.value`, `sig`, `method`.

## Examples

``` r
stat_hov_test(ToothGrowth, .cols = len, group = supp)
#> Variance homogeneity test 
#> # A tibble: 1 × 8
#>   variable n_used statistic    df   df2 p.value method                     sig  
#> * <chr>     <int>     <dbl> <int> <int>   <dbl> <chr>                      <chr>
#> 1 len          60      1.21     1    58   0.275 median-centered Levene te… ""   
stat_hov_test(mtcars, .cols = mpg, group = cyl, method = "bartlett")
#> Variance homogeneity test 
#> # A tibble: 1 × 8
#>   variable n_used statistic    df   df2 p.value method                     sig  
#> * <chr>     <int>     <dbl> <dbl> <dbl>   <dbl> <chr>                      <chr>
#> 1 mpg          32      8.39     2    NA  0.0150 Bartlett test of homogene… *    
```
