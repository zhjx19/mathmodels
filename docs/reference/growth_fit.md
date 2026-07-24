# Nonlinear Growth Curve Fit

Fits nonlinear growth, saturation, and sigmoid curves via
[`minpack.lm::nlsLM()`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html)
with automatic starting value generation.

## Usage

``` r
growth_fit(x, y, type = c("logistic", "gompertz", "saturation", "mm"))
```

## Arguments

- x:

  Numeric vector of observed predictor values.

- y:

  Numeric vector of observed response values.

- type:

  Character, one of `"logistic"`, `"gompertz"`, `"saturation"`, or
  `"mm"`.

## Value

A named list, see
[`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md)
for details. `model` contains the `nls` object.

## Details

Models and starting value strategies:

**Logistic**: `y = L / (1 + exp(-k * (x - x0)))`

- `L0 = max(y) * 1.05`

- `lm(log((L0 - y) / y) ~ x)` for `k, x0` initials

**Gompertz**: `y = L * exp(-exp(-k * (x - x0)))`

- `L0 = max(y) * 1.05`

- `lm(log(-log(y / L0)) ~ x)` for `k, x0` initials

**Saturation** (exponential approach): `y = a * (1 - exp(-b * x))`

- `a0 = max(y) * 1.1`

- `lm(log(a0 - y) ~ x)` for `b` initial

**Michaelis-Menten**: `y = a * x / (b + x)`

- `a0 = max(y) * 1.1`

- Lineweaver-Burk `lm(1/y ~ 1/x)` for `a, b` initials

Falls back to heuristic defaults if linearized estimation fails.

## Examples

``` r
# \donttest{
set.seed(42)
x = seq(0, 10, length.out = 30)
y = 100 / (1 + exp(-1.5 * (x - 5))) + rnorm(30, sd = 2)
growth_fit(x, y, type = "logistic")
#> $model
#> Nonlinear regression model
#>   model: y ~ L/(1 + exp(-k * (x - x0)))
#>    data: parent.frame()
#>       L       k      x0 
#> 100.288   1.373   5.025 
#>  residual sum-of-squares: 151.4
#> 
#> Number of iterations to convergence: 6 
#> Achieved convergence tolerance: 1.49e-08
#> 
#> $coefficient
#> # A tibble: 3 × 7
#>   term  estimate std.error statistic  p.value conf.low conf.high
#>   <chr>    <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 L       100.      0.878      114.  8.65e-38    98.5     102.  
#> 2 k         1.37    0.0555      24.8 4.29e-20     1.26      1.49
#> 3 x0        5.02    0.0347     145.  1.45e-40     4.95      5.10
#> 
#> $model_info
#> # A tibble: 1 × 5
#>   r.squared adj.r.squared   aic   bic sigma
#>       <dbl>         <dbl> <dbl> <dbl> <dbl>
#> 1     0.997         0.997  140.  144.  2.37
#> 
#> $residuals
#> # A tibble: 30 × 3
#>    .observed .fitted .residual
#>        <dbl>   <dbl>     <dbl>
#>  1     2.80    0.101     2.70 
#>  2    -1.04    0.162    -1.20 
#>  3     0.882   0.260     0.622
#>  4     1.53    0.417     1.11 
#>  5     1.24    0.667     0.577
#>  6     0.517   1.07     -0.550
#>  7     4.24    1.70      2.54 
#>  8     1.84    2.71     -0.870
#>  9     7.39    4.27      3.11 
#> 10     5.37    6.69     -1.32 
#> # ℹ 20 more rows
#> 
#> $formula
#> [1] "y = 100.3 / (1 + exp(-1.373 * (x - 5.025)))"
#> 
#> $type
#> [1] "logistic"
#> 
#> $input
#> # A tibble: 30 × 2
#>        x      y
#>    <dbl>  <dbl>
#>  1 0      2.80 
#>  2 0.345 -1.04 
#>  3 0.690  0.882
#>  4 1.03   1.53 
#>  5 1.38   1.24 
#>  6 1.72   0.517
#>  7 2.07   4.24 
#>  8 2.41   1.84 
#>  9 2.76   7.39 
#> 10 3.10   5.37 
#> # ℹ 20 more rows
#> 
# }
```
