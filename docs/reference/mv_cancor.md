# Canonical Correlation Analysis

Relates two blocks of numeric columns via
[`stats::cancor`](https://rdrr.io/r/stats/cancor.html) on standardized
data, reporting canonical correlations with Wilks' lambda significance
tests, raw coefficients, and structure correlations.

## Usage

``` r
mv_cancor(data, x_cols, y_cols, scale = TRUE)
```

## Arguments

- data:

  A data frame.

- x_cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns of the first set.

- y_cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns of the second set; must not overlap `x_cols`.

- scale:

  Logical. If `TRUE` (default), both sets are standardized so
  coefficients live on the standardized scale.

## Value

An object of class `mv_cancor`: a list with components `correlations`
(tibble with `dimension`, `canonical_correlation`, `wilks_lambda`,
`chi_square`, `df`, `p_value`), `x_coef` (tibble with `variable` and
`Can1`...), `y_coef`, `x_structure` (tibble with `variable` and
`Can1`... structure correlations), `y_structure`, and `meta` (where
`n_dropped` counts rows incomplete in *either* column set).

## Examples

``` r
r = mv_cancor(mtcars, mpg:disp, hp:wt)
r$correlations
#> # A tibble: 3 × 6
#>   dimension canonical_correlation wilks_lambda chi_square    df  p_value
#>       <int>                 <dbl>        <dbl>      <dbl> <dbl>    <dbl>
#> 1         1                0.965        0.0514    81.6        9 7.64e-14
#> 2         2                0.510        0.739      8.31       4 8.09e- 2
#> 3         3                0.0320       0.999      0.0281     1 8.67e- 1
```
