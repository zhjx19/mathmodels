# Analysis of Variance (ANOVA)

Fits an ANOVA on a tidy data frame: `group` columns become factors and
are entered with full interactions when several are supplied. One row
per model term (plus a `Residuals` row); `resp` and `.by` are repeated
(row-per-response and row-per-slice).

## Usage

``` r
stat_anova(data, resp, group, .by = NULL)
```

## Arguments

- data:

  A data frame.

- resp:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric response column(s); each column is analysed separately.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Factor column(s); several columns are entered with interactions
  (`g1 * g2`).

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns.

## Value

A tibble with class `stat_infer`: `resp`, `.by` identifiers, `term`,
`df`, `sumsq`, `meansq`, `statistic`, `p.value`, `sig`, `effect`
(partial eta-squared), `method`.

## Details

The table is the sequential decomposition from
[`stats::anova()`](https://rdrr.io/r/stats/anova.html) — Type I
(sequential) sums of squares. With unbalanced designs the early terms
therefore differ from the Type II/III decompositions of
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) or
[`rstatix::anova_test()`](https://rpkgs.datanovia.com/rstatix/reference/anova_test.html)
(see
[`stat_ancova()`](https://zhjx19.github.io/mathmodels/reference/stat_ancova.md)
for a Type II ANCOVA table).

## Examples

``` r
stat_anova(ToothGrowth, resp = len, group = supp)
#> Analysis of variance (ANOVA) 
#> # A tibble: 2 × 10
#>   resp  term         df sumsq meansq statistic p.value  effect p.adjusted sig  
#> * <chr> <chr>     <int> <dbl>  <dbl>     <dbl>   <dbl>   <dbl>      <dbl> <chr>
#> 1 len   supp          1  205.  205.       3.67  0.0604  0.0595     0.0604 "."  
#> 2 len   Residuals    58 3247.   56.0     NA    NA      NA         NA      ""   
stat_anova(mtcars, resp = mpg, group = c(am, cyl))
#> Analysis of variance (ANOVA) 
#> # A tibble: 4 × 10
#>   resp  term         df sumsq meansq statistic  p.value  effect p.adjusted sig  
#> * <chr> <chr>     <int> <dbl>  <dbl>     <dbl>    <dbl>   <dbl>      <dbl> <chr>
#> 1 mpg   am            1 405.  405.       44.1   4.85e-7  0.629     4.85e-7 "***"
#> 2 mpg   cyl           2 456.  228.       24.8   9.35e-7  0.656     9.35e-7 "***"
#> 3 mpg   am:cyl        2  25.4  12.7       1.38  2.69e-1  0.0962    2.69e-1 ""   
#> 4 mpg   Residuals    26 239.    9.19     NA    NA       NA        NA       ""   
```
