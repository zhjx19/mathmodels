# Polynomial Fit

Fits a polynomial regression model via
[`lm()`](https://rdrr.io/r/stats/lm.html).

## Usage

``` r
poly_fit(x, y, degree = 1)
```

## Arguments

- x:

  Numeric vector of observed predictor values.

- y:

  Numeric vector of observed response values.

- degree:

  Integer, degree of the polynomial. Default is 1 (linear).

## Value

A named list with elements: `model` (lm object), `coefficient` (tibble
with estimates, SE, t-value, p-value, 95% CI), `model_info` (tibble with
r.squared, adj.r.squared, AIC, BIC, sigma, all computed on the original
(data) scale), `residuals` (tibble of .observed, .fitted, .residual),
`formula` (character), `type` (character), and `input` (tibble of x, y).

## Examples

``` r
x = 1:20
y = 3 + 2*x + rnorm(20, sd = 2)
poly_fit(x, y, degree = 1)
#> $model
#> 
#> Call:
#> lm(formula = y ~ poly(x, degree, raw = TRUE))
#> 
#> Coefficients:
#>                 (Intercept)  poly(x, degree, raw = TRUE)  
#>                       1.984                        2.075  
#> 
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term                  estimate std.error statistic  p.value conf.low conf.high
#>   <chr>                    <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)               1.98    0.819       2.42 2.63e- 2    0.262      3.71
#> 2 poly(x, degree, raw …     2.07    0.0684     30.3  6.61e-17    1.93       2.22
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.981         0.980  83.4  86.3  1.76
#> 
#> $residuals
#> # A tibble: 20 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1      3.77    4.06   -0.286 
#>  2      2.95    6.13   -3.18  
#>  3      6.55    8.21   -1.66  
#>  4     11.4    10.3     1.08  
#>  5     14.1    12.4     1.78  
#>  6     14.0    14.4    -0.417 
#>  7     17.0    16.5     0.494 
#>  8     21.2    18.6     2.67  
#>  9     23.9    20.7     3.22  
#> 10     20.8    22.7    -1.92  
#> 11     24.8    24.8    -0.0390
#> 12     29.4    26.9     2.52  
#> 13     28.1    29.0    -0.893 
#> 14     30.9    31.0    -0.133 
#> 15     32.8    33.1    -0.275 
#> 16     33.2    35.2    -1.95  
#> 17     36.1    37.3    -1.14  
#> 18     38.9    39.3    -0.385 
#> 19     40.2    41.4    -1.23  
#> 20     45.2    43.5     1.75  
#> 
#> $formula
#> [1] "y = 1.984 + 2.075 * x"
#> 
#> $type
#> [1] "poly"
#> 
#> $input
#> # A tibble: 20 × 2
#>        x     y
#>    <int> <dbl>
#>  1     1  3.77
#>  2     2  2.95
#>  3     3  6.55
#>  4     4 11.4 
#>  5     5 14.1 
#>  6     6 14.0 
#>  7     7 17.0 
#>  8     8 21.2 
#>  9     9 23.9 
#> 10    10 20.8 
#> 11    11 24.8 
#> 12    12 29.4 
#> 13    13 28.1 
#> 14    14 30.9 
#> 15    15 32.8 
#> 16    16 33.2 
#> 17    17 36.1 
#> 18    18 38.9 
#> 19    19 40.2 
#> 20    20 45.2 
#> 
poly_fit(x, y, degree = 2)
#> $model
#> 
#> Call:
#> lm(formula = y ~ poly(x, degree, raw = TRUE))
#> 
#> Coefficients:
#>                  (Intercept)  poly(x, degree, raw = TRUE)1  
#>                      0.68823                       2.42795  
#> poly(x, degree, raw = TRUE)2  
#>                     -0.01683  
#> 
#> 
#> $coefficient
#> # A tibble: 3 × 7
#>   term                   estimate std.error statistic p.value conf.low conf.high
#>   <chr>                     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 (Intercept)              0.688     1.29       0.534 6.00e-1  -2.03      3.41  
#> 2 poly(x, degree, raw =…   2.43      0.283      8.59  1.37e-7   1.83      3.02  
#> 3 poly(x, degree, raw =…  -0.0168    0.0131    -1.29  2.15e-1  -0.0444    0.0108
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.982         0.980  83.5  87.5  1.73
#> 
#> $residuals
#> # A tibble: 20 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1      3.77    3.10   0.673  
#>  2      2.95    5.48  -2.53   
#>  3      6.55    7.82  -1.27   
#>  4     11.4    10.1    1.23   
#>  5     14.1    12.4    1.73   
#>  6     14.0    14.7   -0.636  
#>  7     17.0    16.9    0.141  
#>  8     21.2    19.0    2.21   
#>  9     23.9    21.2    2.70   
#> 10     20.8    23.3   -2.48   
#> 11     24.8    25.4   -0.594  
#> 12     29.4    27.4    2.00   
#> 13     28.1    29.4   -1.35   
#> 14     30.9    31.4   -0.486  
#> 15     32.8    33.3   -0.494  
#> 16     33.2    35.2   -2.00   
#> 17     36.1    37.1   -0.990  
#> 18     38.9    38.9    0.00169
#> 19     40.2    40.7   -0.573  
#> 20     45.2    42.5    2.71   
#> 
#> $formula
#> [1] "y = 0.6882 + 2.428 * x - 0.01683 * x^2"
#> 
#> $type
#> [1] "poly"
#> 
#> $input
#> # A tibble: 20 × 2
#>        x     y
#>    <int> <dbl>
#>  1     1  3.77
#>  2     2  2.95
#>  3     3  6.55
#>  4     4 11.4 
#>  5     5 14.1 
#>  6     6 14.0 
#>  7     7 17.0 
#>  8     8 21.2 
#>  9     9 23.9 
#> 10    10 20.8 
#> 11    11 24.8 
#> 12    12 29.4 
#> 13    13 28.1 
#> 14    14 30.9 
#> 15    15 32.8 
#> 16    16 33.2 
#> 17    17 36.1 
#> 18    18 38.9 
#> 19    19 40.2 
#> 20    20 45.2 
#> 
```
