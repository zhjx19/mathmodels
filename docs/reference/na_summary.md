# Missing Value Summary

Summarizes missing values per column of a tidy data frame: the number
and percentage of missing entries. The row dimension is covered by
[`plot_na_heatmap()`](https://zhjx19.github.io/mathmodels/reference/plot_na_heatmap.md).

## Usage

``` r
na_summary(data, .cols = dplyr::everything())
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to summarize; defaults to all columns.

## Value

A tibble: `variable`, `n_missing`, `pct_missing` (original column
order).

## Examples

``` r
d = data.frame(x = c(1, NA, 3), y = c(NA, NA, 2))
na_summary(d)
#> # A tibble: 2 × 3
#>   variable n_missing pct_missing
#>   <chr>        <int>       <dbl>
#> 1 x                1       0.333
#> 2 y                2       0.667
```
