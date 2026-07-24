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
r.squared, adj.r.squared, AIC, BIC, sigma), `residuals` (tibble of
.observed, .fitted, .residual), `formula` (character), `type`
(character), and `input` (tibble of x, y).

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
#>                       2.535                        2.008  
#> 
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term                  estimate std.error statistic  p.value conf.low conf.high
#>   <chr>                    <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)               2.54    0.941       2.69 1.48e- 2    0.558      4.51
#> 2 poly(x, degree, raw …     2.01    0.0786     25.6  1.35e-15    1.84       2.17
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.973         0.972  88.9  91.9  2.03
#> 
#> $residuals
#> # A tibble: 20 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1      5.91    4.54     1.37 
#>  2      8.41    6.55     1.86 
#>  3     11.1     8.56     2.51 
#>  4      9.78   10.6     -0.784
#>  5     14.0    12.6      1.44 
#>  6     11.6    14.6     -3.02 
#>  7     15.4    16.6     -1.16 
#>  8     17.3    18.6     -1.30 
#>  9     16.2    20.6     -4.43 
#> 10     23.1    22.6      0.460
#> 11     25.4    24.6      0.792
#> 12     26.3    26.6     -0.350
#> 13     30.5    28.6      1.88 
#> 14     29.5    30.6     -1.10 
#> 15     30.3    32.7     -2.39 
#> 16     35.9    34.7      1.21 
#> 17     35.4    36.7     -1.29 
#> 18     41.9    38.7      3.21 
#> 19     40.1    40.7     -0.544
#> 20     44.3    42.7      1.62 
#> 
#> $formula
#> [1] "y = 2.535 + 2.008 * x"
#> 
#> $type
#> [1] "poly"
#> 
#> $input
#> # A tibble: 20 × 2
#>        x     y
#>    <int> <dbl>
#>  1     1  5.91
#>  2     2  8.41
#>  3     3 11.1 
#>  4     4  9.78
#>  5     5 14.0 
#>  6     6 11.6 
#>  7     7 15.4 
#>  8     8 17.3 
#>  9     9 16.2 
#> 10    10 23.1 
#> 11    11 25.4 
#> 12    12 26.3 
#> 13    13 30.5 
#> 14    14 29.5 
#> 15    15 30.3 
#> 16    16 35.9 
#> 17    17 35.4 
#> 18    18 41.9 
#> 19    19 40.1 
#> 20    20 44.3 
#> 
poly_fit(x, y, degree = 2)
#> $model
#> 
#> Call:
#> lm(formula = y ~ poly(x, degree, raw = TRUE))
#> 
#> Coefficients:
#>                  (Intercept)  poly(x, degree, raw = TRUE)1  
#>                      4.88435                       1.36699  
#> poly(x, degree, raw = TRUE)2  
#>                      0.03051  
#> 
#> 
#> $coefficient
#> # A tibble: 3 × 7
#>   term                   estimate std.error statistic p.value conf.low conf.high
#>   <chr>                     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 (Intercept)              4.88      1.37        3.57 2.36e-3  2.00       7.77  
#> 2 poly(x, degree, raw =…   1.37      0.300       4.55 2.81e-4  0.734      2.00  
#> 3 poly(x, degree, raw =…   0.0305    0.0139      2.20 4.21e-2  0.00122    0.0598
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.979         0.977  85.9  89.9  1.84
#> 
#> $residuals
#> # A tibble: 20 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1      5.91    6.28    -0.371
#>  2      8.41    7.74     0.669
#>  3     11.1     9.26     1.81 
#>  4      9.78   10.8     -1.06 
#>  5     14.0    12.5      1.53 
#>  6     11.6    14.2     -2.62 
#>  7     15.4    15.9     -0.517
#>  8     17.3    17.8     -0.475
#>  9     16.2    19.7     -3.49 
#> 10     23.1    21.6      1.47 
#> 11     25.4    23.6      1.80 
#> 12     26.3    25.7      0.596
#> 13     30.5    27.8      2.70 
#> 14     29.5    30.0     -0.456
#> 15     30.3    32.3     -1.99 
#> 16     35.9    34.6      1.30 
#> 17     35.4    36.9     -1.56 
#> 18     41.9    39.4      2.51 
#> 19     40.1    41.9     -1.73 
#> 20     44.3    44.4     -0.117
#> 
#> $formula
#> [1] "y = 4.884 + 1.367 * x + 0.03051 * x^2"
#> 
#> $type
#> [1] "poly"
#> 
#> $input
#> # A tibble: 20 × 2
#>        x     y
#>    <int> <dbl>
#>  1     1  5.91
#>  2     2  8.41
#>  3     3 11.1 
#>  4     4  9.78
#>  5     5 14.0 
#>  6     6 11.6 
#>  7     7 15.4 
#>  8     8 17.3 
#>  9     9 16.2 
#> 10    10 23.1 
#> 11    11 25.4 
#> 12    12 26.3 
#> 13    13 30.5 
#> 14    14 29.5 
#> 15    15 30.3 
#> 16    16 35.9 
#> 17    17 35.4 
#> 18    18 41.9 
#> 19    19 40.1 
#> 20    20 44.3 
#> 
```
