# Logistic Regression

Fits a binary logistic regression model via `glm(family = binomial())`,
with optional stepwise selection.

## Usage

``` r
reg_logistic(
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

- odds_ratio:

  Tibble with odds ratios and 95% CIs.

- residuals:

  Tibble of observed (0/1), fitted probability, Pearson residual.

## Details

Same stepwise logic as
[`reg_lm`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md), but
uses GLM with the binomial family. A Hosmer-Lemeshow goodness-of-fit
test is included in diagnostics.

## Examples

``` r
set.seed(42)
x1 = rnorm(100)
logit = 0.5 + x1
prob = 1 / (1 + exp(-logit))
y = rbinom(100, 1, prob)
reg_logistic(y ~ x1, data = data.frame(y, x1))
#> $model
#> 
#> Call:  glm(formula = formula, family = binomial(), data = data)
#> 
#> Coefficients:
#> (Intercept)           x1  
#>      0.2855       1.0744  
#> 
#> Degrees of Freedom: 99 Total (i.e. Null);  98 Residual
#> Null Deviance:       136.7 
#> Residual Deviance: 114.9     AIC: 118.9
#> 
#> $coefficient
#> # A tibble: 2 × 7
#>   term        estimate std.error statistic   p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)    0.285     0.226      1.26 0.207       -0.163     0.734
#> 2 x1             1.07      0.270      3.98 0.0000682    0.539     1.61 
#> 
#> $odds_ratio
#> # A tibble: 2 × 4
#>   term        odds_ratio conf.low conf.high
#>   <chr>            <dbl>    <dbl>     <dbl>
#> 1 (Intercept)       1.33    0.850      2.08
#> 2 x1                2.93    1.71       5.00
#> 
#> $residuals
#> # A tibble: 100 × 3
#>    observed fitted_val pearson_residual
#>       <int>      <dbl>            <dbl>
#>  1        0      0.853           -2.41 
#>  2        1      0.420            1.17 
#>  3        0      0.663           -1.40 
#>  4        1      0.724            0.617
#>  5        1      0.673            0.698
#>  6        1      0.543            0.918
#>  7        0      0.871           -2.60 
#>  8        1      0.546            0.912
#>  9        1      0.921            0.293
#> 10        1      0.554            0.897
#> # ℹ 90 more rows
#> 
#> $diagnostics
#> # A tibble: 3 × 2
#>   metric     value
#>   <chr>      <dbl>
#> 1 HLR_chisq  9.79 
#> 2 HLR_df     8    
#> 3 HLR_pvalue 0.280
#> 
#> $model_info
#> # A tibble: 1 × 2
#>     aic   bic
#>   <dbl> <dbl>
#> 1  119.  124.
#> 
#> $formula
#> y ~ x1
#> <environment: 0x00000267a2d11968>
#> 
#> $input
#>     y          x1
#> 1   0  1.37095845
#> 2   1 -0.56469817
#> 3   0  0.36312841
#> 4   1  0.63286260
#> 5   1  0.40426832
#> 6   1 -0.10612452
#> 7   0  1.51152200
#> 8   1 -0.09465904
#> 9   1  2.01842371
#> 10  1 -0.06271410
#> 11  1  1.30486965
#> 12  1  2.28664539
#> 13  0 -1.38886070
#> 14  1 -0.27878877
#> 15  1 -0.13332134
#> 16  0  0.63595040
#> 17  0 -0.28425292
#> 18  0 -2.65645542
#> 19  0 -2.44046693
#> 20  0  1.32011335
#> 21  1 -0.30663859
#> 22  0 -1.78130843
#> 23  1 -0.17191736
#> 24  0  1.21467470
#> 25  1  1.89519346
#> 26  0 -0.43046913
#> 27  1 -0.25726938
#> 28  0 -1.76316309
#> 29  1  0.46009735
#> 30  0 -0.63999488
#> 31  1  0.45545012
#> 32  1  0.70483734
#> 33  1  1.03510352
#> 34  0 -0.60892638
#> 35  0  0.50495512
#> 36  1 -1.71700868
#> 37  0 -0.78445901
#> 38  1 -0.85090759
#> 39  0 -2.41420765
#> 40  1  0.03612261
#> 41  1  0.20599860
#> 42  0 -0.36105730
#> 43  1  0.75816324
#> 44  0 -0.72670483
#> 45  0 -1.36828104
#> 46  1  0.43281803
#> 47  0 -0.81139318
#> 48  1  1.44410126
#> 49  1 -0.43144620
#> 50  0  0.65564788
#> 51  1  0.32192527
#> 52  0 -0.78383894
#> 53  1  1.57572752
#> 54  1  0.64289931
#> 55  1  0.08976065
#> 56  0  0.27655075
#> 57  1  0.67928882
#> 58  1  0.08983289
#> 59  0 -2.99309008
#> 60  1  0.28488295
#> 61  0 -0.36723464
#> 62  1  0.18523056
#> 63  1  0.58182373
#> 64  0  1.39973683
#> 65  0 -0.72729206
#> 66  1  1.30254263
#> 67  0  0.33584812
#> 68  1  1.03850610
#> 69  0  0.92072857
#> 70  0  0.72087816
#> 71  0 -1.04311894
#> 72  0 -0.09018639
#> 73  1  0.62351816
#> 74  0 -0.95352336
#> 75  1 -0.54282881
#> 76  1  0.58099650
#> 77  1  0.76817874
#> 78  1  0.46376759
#> 79  0 -0.88577630
#> 80  0 -1.09978090
#> 81  1  1.51270701
#> 82  1  0.25792144
#> 83  1  0.08844023
#> 84  0 -0.12089654
#> 85  0 -1.19432890
#> 86  1  0.61199690
#> 87  1 -0.21713985
#> 88  1 -0.18275671
#> 89  1  0.93334633
#> 90  1  0.82177311
#> 91  1  1.39211638
#> 92  0 -0.47617392
#> 93  1  0.65034856
#> 94  1  1.39111046
#> 95  0 -1.11078888
#> 96  0 -0.86079259
#> 97  1 -1.13173868
#> 98  0 -1.45921400
#> 99  1  0.07998255
#> 100 1  0.65320434
#> 
```
