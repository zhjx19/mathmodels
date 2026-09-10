# Normality Tests

Runs Shapiro-Wilk tests per column (and per `.by` slice), together with
moment skewness and excess kurtosis. Shapiro-Wilk requires n between 3
and 5000 per group; outside that range the result is `NA`.

## Usage

``` r
stat_normality(data, .cols, .by = NULL, p_adjust = "none")
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns (e.g. groups whose normality is checked
  separately).

- p_adjust:

  P-value adjustment across returned rows.

## Value

A tibble with class `stat_infer`: `variable`, `.by` identifiers,
`n_used`, `statistic`, `p.value`, `p.adjusted`, `sig`, `skew`, `kurt`,
`method`.

## Examples

``` r
stat_normality(mtcars, .cols = mpg:hp)
#> Shapiro-Wilk normality test 
#> # A tibble: 4 × 9
#>   variable n_used statistic    p.value   skew    kurt method    p.adjusted sig  
#> * <chr>     <int>     <dbl>      <dbl>  <dbl>   <dbl> <chr>          <dbl> <chr>
#> 1 mpg          32     0.948 0.123       0.640 -0.201  Shapiro-… 0.123      ""   
#> 2 cyl          32     0.753 0.00000606 -0.183 -1.68   Shapiro-… 0.00000606 "***"
#> 3 disp         32     0.920 0.0208      0.400 -1.09   Shapiro-… 0.0208     "*"  
#> 4 hp           32     0.933 0.0488      0.761  0.0522 Shapiro-… 0.0488     "*"  
stat_normality(ToothGrowth, .cols = len, .by = supp)
#> Shapiro-Wilk normality test 
#> # A tibble: 2 × 10
#>   supp  variable n_used statistic p.value   skew   kurt method  p.adjusted sig  
#> * <fct> <chr>     <int>     <dbl>   <dbl>  <dbl>  <dbl> <chr>        <dbl> <chr>
#> 1 OJ    len          30     0.918  0.0236 -0.550 -0.892 Shapir…     0.0236 "*"  
#> 2 VC    len          30     0.966  0.428   0.290 -0.782 Shapir…     0.428  ""   
```
