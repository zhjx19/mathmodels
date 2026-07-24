# Linearizable Curve Fit

Fits common empirical curves via variable transformation and
[`lm()`](https://rdrr.io/r/stats/lm.html). Supported types: exponential,
power-law, logarithmic, and hyperbolic.

## Usage

``` r
curve_fit(x, y, type = c("exp", "power", "log", "hyperbolic"))
```

## Arguments

- x:

  Numeric vector of observed predictor values.

- y:

  Numeric vector of observed response values.

- type:

  Character, one of `"exp"`, `"power"`, `"log"`, `"hyperbolic"`.

## Value

A named list, see
[`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md)
for details.

## Details

Each curve is transformed to a linear form:

- `"exp"`: fits `log(y) ~ x`, back-transforms to `y = a * exp(b * x)`

- `"power"`: fits `log(y) ~ log(x)`, back-transforms to `y = a * x^b`

- `"log"`: fits `y ~ log(x)`, gives `y = a + b * log(x)`

- `"hyperbolic"`: fits `y ~ 1/x`, gives `y = a + b / x`

Fit statistics are computed on the original scale.

**Note on log-transform bias**: For `"exp"` and `"power"` types, fitted
values are back-transformed directly (`a * exp(b * x)` or `a * x^b`),
which yields the geometric-mean prediction. As a result, predicted
values on the original scale are systematically biased downward by a
factor of approximately `exp(sigma^2 / 2)`. This does not affect the
regression coefficients (unbiased on the log scale) or the
goodness-of-fit statistics (computed on the original scale).

**Note on coefficients**: For `"exp"` and `"power"` types, the
coefficient table reports estimates on the log-transformed scale (the
parameterization of the underlying
[`lm()`](https://rdrr.io/r/stats/lm.html) model). The formula string
shows the back-transformed parameters on the original scale.

## Examples

``` r
set.seed(42)
x = 1:20
y = 2 * exp(0.15 * x) * rlnorm(20, sdlog = 0.1)
curve_fit(x, y, type = "exp")
#> $model
#> 
#> Call:
#> lm(formula = ly ~ x)
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.7949       0.1421  
#> 
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.795   0.0586       13.6 6.81e-11    0.672     0.918
#> 2 x              0.142   0.00489      29.1 1.40e-16    0.132     0.152
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.935         0.931 -22.2 -19.2  2.89
#> 
#> $residuals
#> # A tibble: 20 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1      2.67    2.55    0.113 
#>  2      2.55    2.94   -0.391 
#>  3      3.25    3.39   -0.139 
#>  4      3.88    3.91   -0.0274
#>  5      4.41    4.51   -0.0982
#>  6      4.87    5.20   -0.328 
#>  7      6.65    5.99    0.659 
#>  8      6.58    6.90   -0.326 
#>  9      9.44    7.96    1.48  
#> 10      8.91    9.17   -0.266 
#> 11     11.9    10.6     1.29  
#> 12     15.2    12.2     3.02  
#> 13     12.2    14.1    -1.82  
#> 14     15.9    16.2    -0.313 
#> 15     18.7    18.7     0.0536
#> 16     23.5    21.5     1.97  
#> 17     24.9    24.8     0.0870
#> 18     22.8    28.6    -5.78  
#> 19     27.1    33.0    -5.88  
#> 20     45.8    38.0     7.84  
#> 
#> $formula
#> [1] "y = 2.214 * exp(0.1421 * x)"
#> 
#> $type
#> [1] "exp"
#> 
#> $input
#> # A tibble: 20 × 2
#>        x     y
#>    <int> <dbl>
#>  1     1  2.67
#>  2     2  2.55
#>  3     3  3.25
#>  4     4  3.88
#>  5     5  4.41
#>  6     6  4.87
#>  7     7  6.65
#>  8     8  6.58
#>  9     9  9.44
#> 10    10  8.91
#> 11    11 11.9 
#> 12    12 15.2 
#> 13    13 12.2 
#> 14    14 15.9 
#> 15    15 18.7 
#> 16    16 23.5 
#> 17    17 24.9 
#> 18    18 22.8 
#> 19    19 27.1 
#> 20    20 45.8 
#> 
```
