# Defuzzification Methods for Fuzzy Comprehensive Evaluation

Implements defuzzification methods for fuzzy evaluation vectors,
including weighted average and maximum membership methods.

## Usage

``` r
defuzzify(mu, scores, method = "weighted_average")
```

## Arguments

- mu:

  Numeric vector, membership degrees for evaluation levels, in `[0, 1]`.

- scores:

  Numeric vector, scores corresponding to each evaluation level (e.g.,
  c(100, 80, 60, 40) for "Excellent", "Good", "Fair", "Poor").

- method:

  Character, defuzzification method: "weighted_average",
  "max_membership", "centroid".

## Value

Numeric, defuzzified output value.

## Examples

``` r
# Example: Defuzzify fuzzy evaluation vectors for three schemes
mu = c(0.318, 0.351, 0.203, 0.128)
scores = c(30, 60, 75, 90)  # Scores for "Poor", "Fair", "Good", "Excellent"
defuzzify(mu, scores, method = "weighted_average")
#> [1] 57.345
defuzzify(mu, scores, method = "max_membership")
#> [1] 60
defuzzify(mu, scores, method = "centroid")
#> [1] 57.345
```
