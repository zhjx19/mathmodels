# Rank Sum Ratio (RSR) Evaluation

Performs Rank Sum Ratio (RSR) evaluation on a dataset of positive
indicators, computing ranks, weighted RSR values, and a linear
regression model to fit RSR against probit-transformed ranks. Supports
integer or non-integer ranking methods.

## Usage

``` r
rank_sum_ratio(data, w = NULL, method = "int")
```

## Arguments

- data:

  Data frame with positive indicator data; first column is an ID column
  for identifying evaluation objects.

- w:

  Numeric vector, weights for indicators (default = equal weights).

- method:

  Character scalar, ranking method: "int" for integer ranks or "non-int"
  for scaled ranks in `[1, n]` (default = "int").

## Value

A list containing:

- `resultTable`: Data frame with RSR values, ranks, cumulative
  frequencies, probit values, and fitted RSR values.

- `reg`: Linear model object fitting RSR against probit values.

- `rankTable`: Data frame with ranked indicator values.

## Details

The `rank_sum_ratio` function implements the RSR method for evaluating
objects based on positive indicators. It ranks the indicators (using
integer or non-integer methods), computes weighted RSR values, adjusts
ranks with probit transformation, and fits a linear regression model to
relate RSR to probit values. The function assumes the first column of
`data` is an ID column, and weights (`w`) can be provided or set to
equal weights by default.

## Examples

``` r
# Example data
data = data.frame(ID = c("A", "B", "C"), X1 = c(10, 20, 15), X2 = c(5, 10, 8))
w = c(0.4, 0.6)
rank_sum_ratio(data, w, method = "int")
#> $resultTable
#> # A tibble: 3 × 8
#>   ID      RSR  barR     f  sumf barRn Probit RSRfit
#>   <chr> <dbl> <dbl> <int> <int> <dbl>  <dbl>  <dbl>
#> 1 A     0.333     1     1     1 0.333   4.57  0.339
#> 2 C     0.667     2     1     2 0.667   5.43  0.656
#> 3 B     1         3     1     3 0.917   6.38  1.01 
#> 
#> $reg
#> 
#> Call:
#> lm(formula = RSR ~ Probit, data = rltTable)
#> 
#> Coefficients:
#> (Intercept)       Probit  
#>     -1.3389       0.3673  
#> 
#> 
#> $rankTable
#>   ID X1 X2
#> 1  A  1  1
#> 2  B  3  3
#> 3  C  2  2
#> 
```
