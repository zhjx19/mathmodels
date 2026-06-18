# Fuzzy Comprehensive Evaluation

Performs fuzzy comprehensive evaluation using different fuzzy
composition operators to combine factor weights with a fuzzy evaluation
matrix. Suitable for multi-criteria decision analysis with weights from
methods like AHP, entropy, CRITIC, CV, or PCA.

## Usage

``` r
fuzzy_eval(w, R, type)
```

## Arguments

- w:

  Numeric vector, factor weights (e.g., from `combine_weights_linear`).

- R:

  Numeric matrix, fuzzy evaluation matrix with columns as factors and
  rows as evaluation grades. Values should be in `[0, 1]`.

- type:

  Integer or character (1-5), specifying the fuzzy composition operator:

  - 1: Min-max (main factor decisive).

  - 2: Product-max (main factor prominent).

  - 3: Weighted sum (additive average).

  - 4: Bounded sum of mins (min-sum bounded).

  - 5: Normalized min-sum (balanced average).

## Value

A numeric vector of normalized comprehensive evaluation results, summing
to 1.

## Details

The function computes a fuzzy comprehensive evaluation vector `B` based
on the weight vector `w` and fuzzy evaluation matrix `R`. Five
composition operators are supported:

- Type 1 (min-max): `max(min(w, R[j,]))`, emphasizes the main factor.

- Type 2 (product-max): `max(w * R[j,])`, highlights the main factor.

- Type 3 (weighted sum): `sum(w * R[j,])`, additive average.

- Type 4 (bounded sum): `min(1, sum(min(w, R[j,])))`, bounds the sum of
  mins.

- Type 5 (normalized min-sum): `sum(min(w, R[j,]/sum(R[j,])))`, balanced
  average.

The output `B` is normalized to sum to 1. If the sum is zero, an error
is thrown. Uses base R for lightweight implementation.

## Examples

``` r
w = c(0.3, 0.3, 0.3, 0.1)  # weights (e.g., from AHP or entropy)

# fuzzy evaluation matrix (3 grades for 4 factors)
R = matrix(c(0.8, 0.7, 0.6, 0.7,
             0.1, 0.2, 0.2, 0.1,
             0.1, 0.1, 0.2, 0.2), nrow = 3, byrow = TRUE)
# Apply fuzzy comprehensive evaluation
fuzzy_eval(w, R, type = 3)  # Weighted sum
#> [1] 0.70 0.16 0.14
```
