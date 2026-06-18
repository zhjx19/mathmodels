# PCA-Based Weighting Method

Computes indicator weights using Principal Component Analysis (PCA). The
method extracts principal components and uses their variance
contribution to derive objective weights for indicators. Optionally
handles positive/negative directions of indicators, and supports
pre-standardized data.

## Usage

``` r
pca_weight(X, index = NULL, nfs = NULL, varimax = TRUE, method = "abs")
```

## Arguments

- X:

  A numeric data frame or matrix where rows represent samples and
  columns represent indicators.

- index:

  A character vector indicating the direction of each indicator. Use
  `"+"` for positive indicators (higher is better), `"-"` for negative
  indicators (lower is better), and `NA` for already standardized
  indicators (no standardization will be applied).

               If `index = NULL` (default), all indicators are treated as `NA`,
               meaning no standardization is performed.

- nfs:

  Number of principal components to use; by default, all are used.

- varimax:

  Whether to perform Varimax rotation, default is TRUE.

- method:

  Weighting Method. "abs" uses absolute loading values `|a_{ji}|`
  (default), "squared" uses `a_{ji}^2`.

## Value

A list containing:

- w:

  Numeric vector of normalized weights for each indicator.

- s:

  Numeric vector of scores for each sample.

- lambda:

  Eigenvalues of principal components (explained variance).

- varP:

  Proportion of variance explained by selected PCs.

## Examples

``` r
# Example: Using PCA to compute indicator weights
ind = c("+","+","-","-")
pca_weight(iris[1:10, 1:4], ind, nfs = 2)
#> $w
#> [1] 0.3653755 0.2576950 0.2288199 0.1481095
#> 
#> $s
#>  [1]  0.603861299 -0.066459291  0.062381683 -0.570617445  0.562358748
#>  [6]  0.304619966 -0.294827855  0.182708247 -0.777385321 -0.006640031
#> 
#> $cum_contrib
#> [1] 0.8454999
#> 
#> $loading
#>                     PC1        PC2
#> Sepal.Length -0.7089283  0.3581971
#> Sepal.Width  -0.5128833 -0.1963096
#> Petal.Length  0.4706905  0.1075345
#> Petal.Width   0.1132333  0.9064181
#> 
#> $lambda
#> [1] 2.7523578 0.6296417
#> 
#> $varP
#> [1] 0.6880895 0.1574104
#> 
```
