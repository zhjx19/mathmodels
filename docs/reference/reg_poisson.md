# Poisson Regression

Fits a count-response regression model via `glm(family = poisson())`,
with optional stepwise selection.

## Usage

``` r
reg_poisson(
  formula,
  data,
  step = FALSE,
  direction = c("both", "forward", "backward"),
  ...
)
```

## Arguments

- formula:

  A formula, e.g. `y ~ x1 + x2`.

- data:

  A data frame containing the variables in `formula`.

- step:

  Logical. Whether to perform stepwise selection (`TRUE`/`FALSE`).
  Default is `FALSE`.

- direction:

  Character: `"both"` (default), `"forward"`, or `"backward"`.

- ...:

  Additional arguments passed to
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`step()`](https://rdrr.io/r/stats/step.html).

## Value

A named list analogous to
[`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
plus:

- dispersion:

  Dispersion ratio (should be near 1).

- residuals:

  Tibble with observed, fitted, Pearson residual.

## Details

Performs an overdispersion check: dispersion \> 1.5 suggests NB may be
preferable.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
mu = exp(1 + 0.3*x1)
y = rpois(100, mu)
reg_poisson(y ~ x1, data = data.frame(y, x1))
#> $model
#> 
#> Call:  glm(formula = formula, family = poisson(), data = data)
#> 
#> Coefficients:
#> (Intercept)           x1  
#>      0.8483       0.4226  
#> 
#> Degrees of Freedom: 99 Total (i.e. Null);  98 Residual
#> Null Deviance:       167.2 
#> Residual Deviance: 123.8     AIC: 374.2
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.848    0.0687     12.3  5.32e-35    0.712     0.985
#> 2 x1             0.423    0.0664      6.37 1.92e-10    0.291     0.554
#> 
#> $dispersion
#> [1] 1.154
#> 
#> $residuals
#> # A tibble: 100 × 3
#>    observed fitted_val pearson_residual
#>       <int>      <dbl>            <dbl>
#>  1        7       4.17           1.39  
#>  2        2       1.84           0.118 
#>  3        5       2.72           1.38  
#>  4        3       3.05          -0.0297
#>  5        1       2.77          -1.06  
#>  6        2       2.23          -0.156 
#>  7        8       4.42           1.70  
#>  8        2       2.24          -0.163 
#>  9        3       5.48          -1.06  
#> 10        2       2.27          -0.182 
#> # ℹ 90 more rows
#> 
#> $diagnostics
#> # A tibble: 3 × 2
#>   metric      value
#>   <chr>       <dbl>
#> 1 dispersion   1.15
#> 2 aic        374.  
#> 3 bic        379.  
#> 
#> $model_info
#> # A tibble: 1 × 2
#>     aic   bic
#>   <dbl> <dbl>
#> 1  374.  379.
#> 
#> $formula
#> y ~ x1
#> <environment: 0x0000020f090869e0>
#> 
#> $input
#>     y          x1
#> 1   7  1.37095845
#> 2   2 -0.56469817
#> 3   5  0.36312841
#> 4   3  0.63286260
#> 5   1  0.40426832
#> 6   2 -0.10612452
#> 7   8  1.51152200
#> 8   2 -0.09465904
#> 9   3  2.01842371
#> 10  2 -0.06271410
#> 11  4  1.30486965
#> 12  6  2.28664539
#> 13  1 -1.38886070
#> 14  0 -0.27878877
#> 15  2 -0.13332134
#> 16  5  0.63595040
#> 17  3 -0.28425292
#> 18  1 -2.65645542
#> 19  1 -2.44046693
#> 20  7  1.32011335
#> 21  2 -0.30663859
#> 22  0 -1.78130843
#> 23  3 -0.17191736
#> 24  8  1.21467470
#> 25  4  1.89519346
#> 26  4 -0.43046913
#> 27  2 -0.25726938
#> 28  0 -1.76316309
#> 29  1  0.46009735
#> 30  1 -0.63999488
#> 31  2  0.45545012
#> 32  2  0.70483734
#> 33  3  1.03510352
#> 34  0 -0.60892638
#> 35  9  0.50495512
#> 36  3 -1.71700868
#> 37  0 -0.78445901
#> 38  4 -0.85090759
#> 39  1 -2.41420765
#> 40  2  0.03612261
#> 41  1  0.20599860
#> 42  3 -0.36105730
#> 43  1  0.75816324
#> 44  1 -0.72670483
#> 45  2 -1.36828104
#> 46  3  0.43281803
#> 47  1 -0.81139318
#> 48  2  1.44410126
#> 49  2 -0.43144620
#> 50  6  0.65564788
#> 51  2  0.32192527
#> 52  1 -0.78383894
#> 53  3  1.57572752
#> 54  3  0.64289931
#> 55  0  0.08976065
#> 56  4  0.27655075
#> 57  1  0.67928882
#> 58  3  0.08983289
#> 59  1 -2.99309008
#> 60  3  0.28488295
#> 61  3 -0.36723464
#> 62  1  0.18523056
#> 63  2  0.58182373
#> 64  8  1.39973683
#> 65  2 -0.72729206
#> 66  2  1.30254263
#> 67  5  0.33584812
#> 68  2  1.03850610
#> 69  6  0.92072857
#> 70  6  0.72087816
#> 71  0 -1.04311894
#> 72  4 -0.09018639
#> 73  3  0.62351816
#> 74  1 -0.95352336
#> 75  4 -0.54282881
#> 76  2  0.58099650
#> 77  2  0.76817874
#> 78  3  0.46376759
#> 79  2 -0.88577630
#> 80  1 -1.09978090
#> 81  4  1.51270701
#> 82  0  0.25792144
#> 83  1  0.08844023
#> 84  7 -0.12089654
#> 85  1 -1.19432890
#> 86  2  0.61199690
#> 87  2 -0.21713985
#> 88  0 -0.18275671
#> 89  3  0.93334633
#> 90  0  0.82177311
#> 91  6  1.39211638
#> 92  3 -0.47617392
#> 93  2  0.65034856
#> 94  3  1.39111046
#> 95  1 -1.11078888
#> 96  0 -0.86079259
#> 97  3 -1.13173868
#> 98  2 -1.45921400
#> 99  1  0.07998255
#> 100 0  0.65320434
#> 
```
