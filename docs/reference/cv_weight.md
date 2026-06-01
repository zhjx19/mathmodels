# Coefficient of Variation Weighting

Computes weights for indicators using the Coefficient of Variation (CV)
method. Weights are derived by normalizing the CV (standard deviation
divided by mean) for each indicator.

## Usage

``` r
cv_weight(X)
```

## Arguments

- X:

  Numeric matrix or data frame with positive indicator data.

## Value

Numeric vector of weights for the indicators, summing to 1.

## Details

The `cv_weight` function calculates weights using the CV method. For
each column in `data`, the CV is computed as the standard deviation
divided by the mean. Weights are obtained by normalizing the CVs to sum
to 1. This lightweight implementation uses base R and assumes all
columns are numeric indicators.

## Examples

``` r
X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
cv_weight(X)
#>       x1       x2 
#> 0.503839 0.496161 
```
