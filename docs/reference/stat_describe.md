# Descriptive Statistics

Computes per-variable descriptive statistics (optionally grouped by
slices of `.by`): counts, mean, sd, median, IQR, t-based confidence
interval for the mean, and moment skewness / excess kurtosis. Statistics
for each variable use only that variable's own non-missing values
(per-variable deletion); rows where other selected variables are missing
do not affect a variable's results. `n_missing` counts the missing
values of the variable itself within the slice.

## Usage

``` r
stat_describe(data, .cols, .by = NULL, conf_level = 0.95)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; statistics are computed separately for each
  slice.

- conf_level:

  Confidence level for the mean CI.

## Value

A tibble with one row per variable (per slice) and class `stat_infer`:
`variable`, `.by` identifiers, `n_used`, `n_missing`, `mean`, `sd`,
`median`, `iqr`, `ci_low`, `ci_high`, `skew`, `kurt`. `n_used` /
`n_missing` are per-variable within-slice counts (other variables'
missing values do not affect them).

## Examples

``` r
stat_describe(mtcars, .cols = mpg:hp)
#> Descriptive statistics 
#> # A tibble: 4 × 11
#>   variable n_used n_missing   mean     sd median    iqr ci_low ci_high   skew
#> * <chr>     <int>     <int>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>
#> 1 mpg          32         0  20.1    6.03   19.2   7.38  17.9    22.3   0.640
#> 2 cyl          32         0   6.19   1.79    6     4      5.54    6.83 -0.183
#> 3 disp         32         0 231.   124.    196.  205.   186.    275.    0.400
#> 4 hp           32         0 147.    68.6   123    83.5  122.    171.    0.761
#> # ℹ 1 more variable: kurt <dbl>
stat_describe(mtcars, .cols = mpg, .by = cyl)
#> Descriptive statistics 
#> # A tibble: 3 × 12
#>     cyl variable n_used n_missing  mean    sd median   iqr ci_low ci_high   skew
#> * <dbl> <chr>     <int>     <int> <dbl> <dbl>  <dbl> <dbl>  <dbl>   <dbl>  <dbl>
#> 1     4 mpg          11         0  26.7  4.51   26    7.6    23.6    29.7  0.299
#> 2     6 mpg           7         0  19.7  1.45   19.7  2.35   18.4    21.1 -0.199
#> 3     8 mpg          14         0  15.1  2.56   15.2  1.85   13.6    16.6 -0.405
#> # ℹ 1 more variable: kurt <dbl>
```
