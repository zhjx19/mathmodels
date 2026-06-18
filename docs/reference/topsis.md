# TOPSIS Method for Multi-Criteria Decision Making

Implements the Technique for Order of Preference by Similarity to Ideal
Solution (TOPSIS) to rank alternatives based on multiple criteria. The
function normalizes the decision matrix using L2 norm, applies weights,
and computes relative closeness to the ideal solution.

## Usage

``` r
topsis(X, w = NULL, index = NULL)
```

## Arguments

- X:

  A numeric matrix or data frame where rows represent alternatives and
  columns represent criteria.

- w:

  A numeric vector of weights for each criterion. Must be non-negative
  and sum to 1. If not provided, equal weights are used.

- index:

  A character vector indicating the direction of each indicator: Use
  `"+"` for positive indicators (higher is better), `"-"` for negative
  indicators (lower is better). If `index = NULL` (default), all
  indicators are treated as `"+"`.

## Value

A named numeric vector of relative closeness scores (in `[0, 1]`) for
each alternative. Higher values indicate better alternatives. Names are
taken from `rownames(X)` or default to "Sample1", "Sample2", etc.

## Details

The TOPSIS method ranks alternatives by:

1.  Normalizing the decision matrix using L2 norm normalization.

2.  Applying weights to form a weighted normalized matrix.

3.  Identifying positive and negative ideal solutions based on indicator
    directions.

4.  Computing Euclidean distances to ideal solutions.

5.  Calculating relative closeness as `S0 / (S0 + Sstar)`, where `S0` is
    the distance to the negative ideal and `Sstar` is the distance to
    the positive ideal.

This implementation supports both positive and negative indicators via
the `index` parameter.

## Examples

``` r
A = data.frame(
  X1 = c(2, 5, 3),  # "+"
  X2 = c(8, 1, 6)   # "-"
)
w = c(0.6, 0.4)
idx = c("+","-")
topsis(A, w, idx)
#> [1] 0.0000000 1.0000000 0.3111391
```
