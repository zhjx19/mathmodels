# Grey Prediction Models

Implements grey prediction models for time series forecasting: `GM11`
applies the GM(1,1) model with level ratio test. `GM1N` applies the
GM(1,N) model with multiple related factors. `DGM21` applies the
DGM(2,1) model for second-order dynamics. `verhulst` applies the
Verhulst model for logistic growth.

## Usage

``` r
GM11(X)

GM1N(dat, new_data = NULL)

DGM21(X)

verhulst(X)
```

## Arguments

- X:

  For `GM11`, `DGM21`, `verhulst`: Numeric vector of original time
  series data.

- dat:

  For `GM1N`: Data frame or matrix, last column is characteristic
  series, others are related factors.

- new_data:

  For `GM1N`: Optional future values of related factors (1 row, m-1
  columns).

## Value

For `GM11`: List with fitted values (`fitted`), next prediction
(`pnext`), prediction function (`f`), matrix (`mat`), parameters (`u`),
level ratios (`lambda`), and range (`rng`). For `GM1N`: List with fitted
values (`fitted`), posterior variance ratio (`C`), small error
probability (`P`), and prediction function (`f`). For `DGM21`,
`verhulst`: List with fitted values (`fitted`), next prediction
(`pnext`), prediction function (`f`), matrix (`mat`), and parameters
(`u`).

## Examples

``` r
# Sample time series for GM11, DGM21, Verhulst
x = c(100, 120, 145, 175, 210)

# GM11
result = GM11(x)
#> Level ratio test passed!
result$fitted    # Fitted values
#> [1] 100.0000 119.9184 144.3784 173.8275 209.2834
result$pnext     # Next prediction
#> [1] 251.9713
result$f(6:8)    # Predict next 3 periods
#> [1] 251.9713 303.3663 365.2444

# DGM21
x = c(2.874,3.278,3.39,3.679,3.77,3.8)
result = DGM21(x)
result$fitted    # Fitted values
#> [1] 2.874000 3.086001 3.408835 3.620105 3.758366 3.848848
result$pnext     # Next prediction
#> [1] 3.908061
result$f(6:8)    # Predict next 3 periods
#> [1] 20.59616 24.50422 28.45103

# Verhulst
x = c(4.93,2.33,3.87,4.35,6.63,7.15,5.37,6.39,7.81,8.35)
result = verhulst(x)
result$fitted    # Fitted values
#>  [1] 4.930000 1.952177 2.635709 3.481640 4.468640 5.528334 6.536384 7.326754
#>  [9] 7.737430 7.673378
result$pnext     # Next prediction
#> [1] 7.149904
result$f(6:8)    # Predict next 3 periods
#> [1] 22.99650 29.53288 36.85964

# Sample data for GM1N
data = data.frame(
  factor1 = c(50, 55, 60, 65, 70),
  factor2 = c(20, 22, 25, 28, 30),
  output = c(100, 120, 145, 175, 210)
)
result = GM1N(data)
#> Smoothness check: Passed 
#> Level ratio check: Passed 
result$fitted
#> [1] 100.0000 111.3703 195.8204 231.7904 301.9012
```
