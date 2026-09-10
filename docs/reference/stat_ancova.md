# Analysis of Covariance (ANCOVA)

Fits `resp ~ group + covariates` (group as factor, covariates numeric)
and returns the ANOVA table with partial eta-squared. The table uses
Type II sums of squares (via
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html)) so the group
term is covariate-adjusted — the ANCOVA question proper. Optionally
checks the parallel-slopes assumption (group-by-covariate interactions)
for every response.

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

  Logical; if `TRUE`, the group-by-covariate interactions are tested for
  every response and a warning (listing all violating responses) is
  issued when any is significant at 0.05 (parallel-slopes assumption in
  doubt).

## Value

A tibble with class `stat_infer`: `resp`, `.by` identifiers, `term`,
`df`, `sumsq`, `meansq`, `statistic`, `p.value`, `sig`, `effect`
(partial eta-squared), `method`. Sums of squares are Type II
(covariate-adjusted group term).

## Details

The model is fitted with
[`stats::aov()`](https://rdrr.io/r/stats/aov.html) and decomposed with
`car::Anova(fit, type = 2)`, so each term is adjusted for all other
terms — for unbalanced group sizes the group sum of squares differs from
the sequential (Type I) decomposition of
[`stats::anova()`](https://rdrr.io/r/stats/anova.html). `test_parallel`
refits `resp ~ group * covariates` per response and covariate and warns
once, listing every response whose interaction is significant.

## Examples

``` r
set.seed(1)
d = data.frame(
  g = rep(c("a", "b"), each = 10),
  x = rep(1:10, 2),
  y = 2 * rep(1:10, 2) + rep(c(0, 5), each = 10) + rnorm(20)
)
stat_ancova(d, resp = y, group = g, covariates = x)
#> Analysis of covariance (ANCOVA) 
#> # A tibble: 3 × 9
#>   resp  term         df sumsq  meansq statistic   p.value effect sig  
#> * <chr> <chr>     <dbl> <dbl>   <dbl>     <dbl>     <dbl>  <dbl> <chr>
#> 1 y     g             1 131.  131.         145.  9.50e-10  0.895 "***"
#> 2 y     x             1 695.  695.         770.  1.35e-15  0.978 "***"
#> 3 y     Residuals    17  15.3   0.902       NA  NA        NA     ""   
```
