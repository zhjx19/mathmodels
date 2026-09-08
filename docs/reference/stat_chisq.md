# Chi-square Test of Independence

Tests the independence of two factor columns. Reports Cramer's V and
warns when expected frequencies are low (suggesting `exact = TRUE`).

## Usage

``` r
stat_chisq(data, row, col, .by = NULL, correct = TRUE, exact = FALSE)
```

## Arguments

- data:

  A data frame.

- row:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Single factor column.

- col:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Single factor column, different from `row`.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

- correct:

  Logical; continuity correction for 2x2 tables.

- exact:

  Logical; Fisher's exact test instead of chi-square.

## Value

A tibble with class `stat_infer`: `.by` identifiers, `n_used`,
`statistic`, `df`, `p.value`, `sig`, `effect` (Cramer's V),
`min_expected`, `pct_below_5`, `method`.

## Examples

``` r
d = data.frame(
  treat = rep(c("a", "b"), each = 30),
  outcome = sample(c("yes", "no"), 60, TRUE)
)
stat_chisq(d, row = treat, col = outcome)
#> Chi-square test of independence 
#> # A tibble: 1 × 9
#>   n_used statistic    df p.value effect min_expected pct_below_5 method    sig  
#> *  <int>     <dbl> <int>   <dbl>  <dbl>        <dbl>       <dbl> <chr>     <chr>
#> 1     60         0     1       1      0           13           0 Pearson'… ""   
```
