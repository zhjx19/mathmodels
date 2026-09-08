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
#>    supp variable n_used statistic    p.value       skew       kurt
#> OJ   OJ      len     30 0.9178434 0.02358773 -0.5504998 -0.8923989
#> VC   VC      len     30 0.9656715 0.42844646  0.2900364 -0.7820183
#>                         method p.adjusted sig
#> OJ Shapiro-Wilk normality test 0.02358773   *
#> VC Shapiro-Wilk normality test 0.42844646    
```
