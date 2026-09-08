# Kruskal-Wallis Rank Sum Test

Nonparametric one-way ANOVA comparing a numeric response across the
levels of a single factor column; effect size is epsilon-squared.

## Usage

``` r
stat_kruskal(data, resp, group, .by = NULL)
```

## Arguments

- data:

  A data frame.

- resp:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric response column(s); each column is analysed separately.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Single factor column (two or more levels).

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

## Value

A tibble with class `stat_infer`: `resp`, `.by` identifiers, `n_used`,
`statistic`, `df`, `p.value`, `effect` (epsilon-squared), `method`.

## Examples

``` r
stat_kruskal(ToothGrowth, resp = len, group = supp)
#> Kruskal-Wallis rank sum test 
#> # A tibble: 1 × 8
#>   resp  n_used statistic    df p.value effect method                       sig  
#> * <chr>  <int>     <dbl> <int>   <dbl>  <dbl> <chr>                        <chr>
#> 1 len       60      3.45     1  0.0634 0.0422 Kruskal-Wallis rank sum test .    
```
