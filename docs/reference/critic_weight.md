# CRITIC Weight Method

Computes objective weights of indicators and scores of samples using the
CRITIC method. The method considers both the contrast intensity (e.g.,
standard deviation or entropy) and conflict among indicators (based on
correlation) to determine indicator importance. This version supports
different methods for contrast intensity and correlation types.

## Usage

``` r
critic_weight(
  X,
  index = NULL,
  method = "std",
  cor_method = "pearson",
  epsilon = 0.002
)
```

## Arguments

- X:

  A numeric data frame or matrix where rows represent samples
  (observations) and columns represent indicators (variables).

- index:

  A character vector indicating the direction of each indicator: Use
  `"+"` for positive indicators (higher is better), `"-"` for negative
  indicators (lower is better), and `NA` for already normalized
  indicators (no rescaling will be applied).

               If `index = NULL` (default), all indicators are treated as `NA`, meaning no normalization is performed.

- method:

  Character scalar; specifies the method used to compute contrast
  intensity. Options: `"std"` (standard deviation, default), or
  `"entropy"` (based on information redundancy).

- cor_method:

  Character scalar; specifies the method for computing correlations.
  Options: `"pearson"` (default), `"spearman"`, or `"kendall"`.

- epsilon:

  A small constant used to replace exact 0s and 1s in the data to
  prevent log(0) errors. Default is 0.002. Only used when
  `method = "entropy"`.

## Value

A list containing:

- w:

  Numeric vector of weights for each indicator.

- s:

  Numeric vector of scores for each sample (row), scaled by 100.

## Examples

``` r
# Example: Using CRITIC method on a simple dataset
X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
critic_weight(X, index)
#> $w
#>        x1        x2 
#> 0.5075187 0.4924813 
#> 
#> $s
#> [1] 59.39851 46.86716 32.83209 50.75187
#> 
critic_weight(X, index, method = "entropy")
#> $w
#>        x1        x2 
#> 0.5455552 0.4544448 
#> 
#> $s
#> [1] 56.26470 47.88147 30.40543 54.53730
#> 
```
