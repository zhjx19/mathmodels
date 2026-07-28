# Grey-Markov Prediction Model

Fits a GM(1,1) model and then uses a Markov chain on relative error
states to correct both historical fitted values and future predictions.

## Usage

``` r
GM11_markov(X, n_ahead = 3, breaks = c(-0.1, -0.05, -0.02, 0, 0.02, 0.05, 0.1))
```

## Arguments

- X:

  Numeric vector of original time series data.

- n_ahead:

  Number of future periods to predict. Must be a positive integer.

- breaks:

  Numeric vector of boundaries for classifying relative errors into
  states. `-Inf` and `Inf` are prepended/appended internally.

## Value

A data frame with columns:

- Period:

  Time period labels, such as `T1`, `T2`, ..., `T+n`.

- Raw:

  Original values; future periods are `NA`.

- GM11_fitted:

  GM(1,1) fitted or predicted values.

- err_state:

  Relative error state for historical periods, and predicted state for
  future periods.

- adj_eff:

  Adjustment effect, defined as the midpoint of the state interval.

- Markov_adj:

  Markov-chain-adjusted fitted or predicted values.

## Details

The function first fits a GM(1,1) model to the original series,
classifies relative errors into states using `breaks`, builds a Markov
chain on these error states, and finally adjusts the GM(1,1) fitted and
predicted values.

## Examples

``` r
X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
GM11_markov(X, n_ahead = 3)
#> Level ratio test passed!
#>    Period   Raw GM11_fitted err_state adj_eff Markov_adj
#> 1      T1 174.0      174.00         5   0.010     174.00
#> 2      T2 179.0      172.81         6   0.035     174.55
#> 3      T3 183.0      183.94         4  -0.010     182.11
#> 4      T4 189.0      195.78         3  -0.035     189.16
#> 5      T5 207.0      208.38         4  -0.010     206.32
#> 6      T6 234.0      221.80         7   0.075     214.30
#> 7      T7 220.5      236.08         2  -0.075     219.61
#> 8      T8 256.0      251.28         5   0.010     253.82
#> 9      T9 270.0      267.46         5   0.010     270.16
#> 10    T10 285.0      284.68         5   0.010     287.56
#> 11    T+1    NA      303.01         5   0.010     306.07
#> 12    T+2    NA      322.52         5   0.010     325.78
#> 13    T+3    NA      343.29         5   0.010     346.76
```
