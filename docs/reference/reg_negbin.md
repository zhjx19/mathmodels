# Negative Binomial Regression

Fits a negative binomial regression model (for over-dispersed counts)
via MASS::glm.nb().

## Usage

``` r
reg_negbin(formula, data, step = FALSE, ...)
```

## Arguments

- formula:

  A formula, e.g. `count ~ x1 + x2`.

- data:

  A data frame.

- step:

  Logical. Always `FALSE` for negative binomial.

- ...:

  Additional arguments passed to
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Value

A named list analogous to
[`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
plus:

- theta:

  Estimated dispersion parameter for the negative binomial distribution.

- residuals:

  Tibble of observed, fitted, Pearson residual.

## Details

**Note:** The `glm.nb` family does not support stepwise selection. Set
`step = FALSE`.

Unlike
[`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
the negative binomial model handles overdispersion through the estimated
dispersion parameter `theta` (fitted by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)): a small
`theta` indicates strong overdispersion, while `theta` growing large
approaches the Poisson model.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
mu = exp(1 + 0.3*x1)
y = rnbinom(100, size = 1, mu = mu)
reg_negbin(y ~ x1, data = data.frame(y, x1))
#> $model
#> 
#> Call:  MASS::glm.nb(formula = formula, data = data, init.theta = 1.065048329, 
#>     link = log)
#> 
#> Coefficients:
#> (Intercept)           x1  
#>      0.9066       0.3223  
#> 
#> Degrees of Freedom: 99 Total (i.e. Null);  98 Residual
#> Null Deviance:       118 
#> Residual Deviance: 107.8     AIC: 426.8
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.907     0.117      7.74 1.01e-14   0.674      1.14 
#> 2 x1             0.322     0.116      2.78 5.52e- 3   0.0918     0.553
#> 
#> $theta
#> [1] 1.065
#> 
#> $residuals
#> # A tibble: 100 × 3
#>    observed fitted_val pearson_residual
#>       <dbl>      <dbl>            <dbl>
#>  1       10       3.85            1.46 
#>  2        3       2.06            0.380
#>  3        1       2.78           -0.562
#>  4        0       3.04           -0.888
#>  5        1       2.82           -0.568
#>  6        1       2.39           -0.500
#>  7        0       4.03           -0.918
#>  8        2       2.40           -0.144
#>  9       33       4.75            5.55 
#> 10        0       2.43           -0.860
#> # ℹ 90 more rows
#> 
#> $diagnostics
#> # A tibble: 3 × 2
#>   metric  value
#>   <chr>   <dbl>
#> 1 aic    427.  
#> 2 bic    435.  
#> 3 theta    1.06
#> 
#> $model_info
#> # A tibble: 1 × 2
#>     aic   bic
#>   <dbl> <dbl>
#> 1  427.  435.
#> 
#> $formula
#> y ~ x1
#> <environment: 0x000001ecd1e58000>
#> 
#> $input
#>      y          x1
#> 1   10  1.37095845
#> 2    3 -0.56469817
#> 3    1  0.36312841
#> 4    0  0.63286260
#> 5    1  0.40426832
#> 6    1 -0.10612452
#> 7    0  1.51152200
#> 8    2 -0.09465904
#> 9   33  2.01842371
#> 10   0 -0.06271410
#> 11   0  1.30486965
#> 12   7  2.28664539
#> 13   2 -1.38886070
#> 14   0 -0.27878877
#> 15   1 -0.13332134
#> 16   1  0.63595040
#> 17   1 -0.28425292
#> 18   1 -2.65645542
#> 19   1 -2.44046693
#> 20  15  1.32011335
#> 21   0 -0.30663859
#> 22   2 -1.78130843
#> 23   1 -0.17191736
#> 24   3  1.21467470
#> 25   1  1.89519346
#> 26   0 -0.43046913
#> 27   3 -0.25726938
#> 28   3 -1.76316309
#> 29   1  0.46009735
#> 30   7 -0.63999488
#> 31   2  0.45545012
#> 32   1  0.70483734
#> 33   0  1.03510352
#> 34   0 -0.60892638
#> 35   2  0.50495512
#> 36   5 -1.71700868
#> 37   0 -0.78445901
#> 38   2 -0.85090759
#> 39   1 -2.41420765
#> 40   6  0.03612261
#> 41   9  0.20599860
#> 42   3 -0.36105730
#> 43   3  0.75816324
#> 44   0 -0.72670483
#> 45   3 -1.36828104
#> 46   0  0.43281803
#> 47   6 -0.81139318
#> 48   2  1.44410126
#> 49   2 -0.43144620
#> 50   1  0.65564788
#> 51   2  0.32192527
#> 52   1 -0.78383894
#> 53   3  1.57572752
#> 54   1  0.64289931
#> 55   0  0.08976065
#> 56   5  0.27655075
#> 57   0  0.67928882
#> 58   4  0.08983289
#> 59   3 -2.99309008
#> 60   0  0.28488295
#> 61   2 -0.36723464
#> 62   0  0.18523056
#> 63  16  0.58182373
#> 64   7  1.39973683
#> 65   4 -0.72729206
#> 66   1  1.30254263
#> 67   3  0.33584812
#> 68   3  1.03850610
#> 69   0  0.92072857
#> 70   3  0.72087816
#> 71   0 -1.04311894
#> 72   4 -0.09018639
#> 73   7  0.62351816
#> 74   2 -0.95352336
#> 75   2 -0.54282881
#> 76   1  0.58099650
#> 77   0  0.76817874
#> 78   2  0.46376759
#> 79   1 -0.88577630
#> 80   0 -1.09978090
#> 81   2  1.51270701
#> 82   2  0.25792144
#> 83   0  0.08844023
#> 84   2 -0.12089654
#> 85   0 -1.19432890
#> 86   2  0.61199690
#> 87   0 -0.21713985
#> 88   2 -0.18275671
#> 89   4  0.93334633
#> 90   0  0.82177311
#> 91   3  1.39211638
#> 92   4 -0.47617392
#> 93   9  0.65034856
#> 94   1  1.39111046
#> 95   1 -1.11078888
#> 96   2 -0.86079259
#> 97   1 -1.13173868
#> 98   3 -1.45921400
#> 99   0  0.07998255
#> 100  4  0.65320434
#> 
```
