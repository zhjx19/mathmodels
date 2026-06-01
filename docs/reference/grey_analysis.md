# Grey Relational Analysis Functions

A collection of functions for performing grey relational analysis,
including calculation of grey correlation degree and evaluation based on
grey correlation. These functions are designed for decision-making and
data analysis by measuring the relational degree between sequences.

## Usage

``` r
grey_corr(ref, cmp, rho = 0.5, w = NULL)

grey_corr_topsis(X, w, index = NULL, rho = 0.5)
```

## Arguments

- ref:

  Numeric vector, the reference sequence for `grey_corr`.

- cmp:

  Numeric matrix or data frame, the comparison sequences for
  `grey_corr`.

- rho:

  Numeric scalar, the distinguishing coefficient (default = 0.5).

- w:

  Numeric vector, weights for weighted correlation (default = equal
  weights).

- X:

  Numeric matrix or data frame, the decision matrix for
  `grey_corr_topsis`.

- index:

  Character vector indicating indicator direction: Use `"+"` for
  positive indicators (higher is better), `"-"` for negative indicators
  (lower is better), and `NA` for already rescaled indicators (no
  rescaling will be applied). If `index = NULL` (default), all
  indicators are treated as `NA`, meaning no rescaling is performed.

## Value

- grey_corr:

  Returns a numeric vector of grey correlation degrees for each
  comparison sequence.

- grey_corr_topsis:

  Returns a numeric vector of relative closeness (grey correlation
  degrees).

## Details

These functions implement grey relational analysis for evaluating
relationships between sequences or decision alternatives:

- grey_corr:

  Computes the grey correlation degree between a reference sequence
  (`ref`) and comparison sequences (`cmp`) using the distinguishing
  coefficient (`rho`) and optional weights (`w`).

- grey_corr_topsis:

  Evaluates a decision matrix (`X`) by normalizing it, applying weights
  (`w`), computing grey correlation with the ideal sequence. Direction
  of indicators can be specified via `index`.

## Examples

``` r
# Grey correlation degree
ref = c(0.9, 0.8, 0.7)
cmp = data.frame(
  x1 = c(0.9, 0.7, 0.8),
  x2 = c(0.8, 0.9, 0.7),
  x3 = c(0.7, 0.8, 0.9)
)
grey_corr(ref, cmp, rho = 0.5)
#>        x1        x2        x3 
#> 0.6666667 0.6666667 0.5555556 

# Grey correlation evaluation
X = data.frame(x1 = c(8, 7, 6), x2 = c(150, 180, 200), x3 = c(60, 80, 100))
w = c(0.3, 0.4, 0.3)
idx = c("+", "+", "+")
grey_corr_topsis(X, w, idx, rho = 0.5)
#> [1] 0.4642511 0.4421740 0.4642511
```
