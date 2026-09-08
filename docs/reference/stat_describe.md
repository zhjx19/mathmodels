# Descriptive Statistics

Computes per-variable descriptive statistics (optionally grouped by
slices of `.by`): counts, mean, sd, median, IQR, t-based confidence
interval for the mean, and moment skewness / excess kurtosis.

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
`median`, `iqr`, `ci_low`, `ci_high`, `skew`, `kurt`.

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
#>   cyl variable n_used n_missing     mean       sd median  iqr   ci_low  ci_high
#> 4   4      mpg     11        21 26.66364 4.509828   26.0 7.60 23.63389 29.69338
#> 6   6      mpg      7        25 19.74286 1.453567   19.7 2.35 18.39853 21.08718
#> 8   8      mpg     14        18 15.10000 2.560048   15.2 1.85 13.62187 16.57813
#>         skew       kurt
#> 4  0.2990323 -1.3604644
#> 6 -0.1994981 -1.5122667
#> 8 -0.4054773 -0.1765741
```
