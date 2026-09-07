# Partial Correlation Matrix

Computes first-order-controlled partial correlations among the selected
numeric columns by inverting the correlation matrix, with t-tests for
each pairwise partial correlation. Each partial correlation controls for
the other `p - 2` variables and is tested with `df = n - p` degrees of
freedom.

## Usage

``` r
mv_pcor(data, cols)
```

## Arguments

- data:

  A data frame.

- cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns; the set must not be (near-)singular.

## Value

An object of class `mv_pcor`: a list with components `partial_corr`
(long tibble with `variable_1`, `variable_2`, `partial_corr`, `t_stat`,
`df`, `p_value`), `pcor_matrix` (wide tibble with `variable` and one
column per variable), and `meta`.

## Examples

``` r
r = mv_pcor(mtcars, mpg:wt)
r$partial_corr
#> # A tibble: 15 × 6
#>    variable_1 variable_2 partial_corr t_stat    df  p_value
#>    <chr>      <chr>             <dbl>  <dbl> <int>    <dbl>
#>  1 mpg        cyl              -0.290 -1.55     26 0.134   
#>  2 mpg        disp              0.200  1.04     26 0.308   
#>  3 mpg        hp               -0.334 -1.81     26 0.0821  
#>  4 mpg        drat              0.133  0.685    26 0.500   
#>  5 mpg        wt               -0.562 -3.47     26 0.00184 
#>  6 cyl        disp              0.493  2.89     26 0.00766 
#>  7 cyl        hp                0.376  2.07     26 0.0489  
#>  8 cyl        drat             -0.321 -1.73     26 0.0958  
#>  9 cyl        wt               -0.257 -1.35     26 0.187   
#> 10 disp       hp                0.276  1.47     26 0.154   
#> 11 disp       drat             -0.117 -0.601    26 0.553   
#> 12 disp       wt                0.628  4.11     26 0.000348
#> 13 hp         drat              0.402  2.24     26 0.0340  
#> 14 hp         wt               -0.202 -1.05     26 0.302   
#> 15 drat       wt               -0.129 -0.665    26 0.512   
```
