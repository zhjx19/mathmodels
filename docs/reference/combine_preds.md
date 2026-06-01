# Combine Multiple Prediction Results

Combines multiple prediction results (e.g., from grey prediction, time
series, or machine learning models) into a single prediction using a
similarity-based weighting approach, improving prediction accuracy.

## Usage

``` r
combine_preds(x)
```

## Arguments

- x:

  Numeric vector, prediction results to be combined (length \>= 2).

## Value

A list with two elements:

- `a`: Numeric, the combined prediction value.

- `w`: Numeric vector, weights for each prediction in `x`, summing to 1.

## Details

The function combines prediction results by constructing a similarity
matrix based on cosine transformation of pairwise differences. Weights
are derived from the principal eigenvector of the similarity matrix,
ensuring predictions closer to each other have higher influence. For two
predictions, equal weights (0.5, 0.5) are used. If all predictions are
identical, equal weights are assigned. Compatible with the `mathmodels`
package for enhancing prediction models, including grey prediction, time
series, or ensemble machine learning.

## Examples

``` r
# Example: Combine three prediction results
preds = c(100, 102, 98)  # E.g., from grey prediction, ARIMA, or ML models
combine_preds(preds)
#> $a
#> [1] 100
#> 
#> $w
#> [1] 0.4142136 0.2928932 0.2928932
#> 
```
