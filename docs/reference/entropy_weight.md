# Entropy Weight Method

Computes the weights of indicators and scores of samples based on the
entropy method. This method objectively determines the importance of
each indicator according to the amount of information it contains.

## Usage

``` r
entropy_weight(X, index = NULL, epsilon = 0.002)
```

## Arguments

- X:

  A numeric data frame or matrix where rows represent samples
  (observations) and columns represent indicators (variables).

- index:

  A character vector indicating the direction of each indicator. Use
  `"+"` for positive indicators (higher is better), `"-"` for negative
  indicators (lower is better), and `NA` for already normalized
  indicators (no rescaling will be applied, but minor adjustments will
  still be made to avoid log(0) errors). If `index = NULL` (default),
  all indicators are treated as `NA`, meaning no normalization or
  rescaling is performed, but a small adjustment is still applied to
  prevent log(0) errors.

- epsilon:

  A small constant used to replace exact 0s and 1s in the data to
  prevent log(0) errors. Default is 0.002.

## Value

A list containing:

- w:

  Numeric vector of weights for each indicator.

- s:

  Numeric vector of scores for each sample (row), scaled by 100.

## Examples

``` r
X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
entropy_weight(X, index)
#> $w
#>        x1        x2 
#> 0.5452618 0.4547382 
#> 
#> $s
#> [1] 56.35354 47.88215 30.39461 54.50808
#> 
```
