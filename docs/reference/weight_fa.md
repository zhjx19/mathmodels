# Factor-Analysis-Based Weighting Method

Computes indicator weights via exploratory factor analysis
([`stats::factanal`](https://rdrr.io/r/stats/factanal.html)). The
(optionally rotated) loading matrix is combined with the variance share
of each factor to derive objective weights: \$\$w_j \propto \sum_i
\mathrm{share}\_i \cdot \|a\_{ji}\|\$\$ where \\a\_{ji}\\ is the loading
of indicator \\j\\ on factor \\i\\ and `share_i` is the proportion of
total variance explained by factor \\i\\. Optionally handles
positive/negative indicator directions and supports pre-standardized
data.

## Usage

``` r
weight_fa(X, index = NULL, nfs = 1, rotation = "varimax", method = "abs")
```

## Arguments

- X:

  A numeric data frame or matrix where rows represent samples and
  columns represent indicators.

- index:

  A character vector indicating the direction of each indicator. Use
  `"+"` for positive indicators (higher is better), `"-"` for negative
  indicators (lower is better), and `NA` for already standardized
  indicators (no standardization will be applied). If `index = NULL`
  (default), all indicators are treated as `NA`, meaning no
  standardization is performed.

- nfs:

  Number of common factors to extract; must be at least 1 and smaller
  than the number of indicators (default 1).

- rotation:

  Rotation passed to
  [`stats::factanal`](https://rdrr.io/r/stats/factanal.html):
  `"varimax"` (default), `"promax"`, `"none"`, or the name of any
  available rotation function.

- method:

  Weighting method: `"abs"` uses absolute loadings \\\|a\_{ji}\|\\
  (default), `"squared"` uses \\a\_{ji}^2\\.

## Value

A list containing:

- w:

  Numeric vector of normalized weights for each indicator.

- s:

  Numeric vector of scores for each sample.

- loading:

  The (rotated) loading matrix.

- share:

  Variance share of each retained factor.

- uniqueness:

  Uniquenesses of the indicators.

## Examples

``` r
set.seed(1)
X = data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
weight_fa(X, nfs = 1)
#> $w
#> [1] 0.93834914 0.03686116 0.02478970
#> 
#> $s
#>  [1] -0.58853646  0.15080559 -0.79411880  1.45921823  0.34578929 -0.65307545
#>  [7]  0.46160939  0.67688129  0.57080830 -0.24983856  1.49134526  0.35291860
#> [13] -0.52201013 -2.09326012  1.02303915 -0.04494190 -0.08965751  0.93275081
#> [19]  0.78849208  0.63297581  0.86730693  0.74104101  0.08716071 -1.90558953
#> [25]  0.53291888 -0.02425853 -0.16435447 -1.38097144 -0.46282955  0.36240612
#> [31]  1.25544537 -0.11603213  0.42037261 -0.14428888 -1.26266963 -0.41522499
#> [37] -0.33825566 -0.07996515  1.02968195  0.72456043 -0.22183489 -0.16405116
#> [43]  0.65549460  0.53666434 -0.61546412 -0.66190360  0.34678874  0.70045222
#> [49] -0.18244761  0.76866647
#> 
#> $loading
#>        Factor1
#> x1  0.99749687
#> x2 -0.03918466
#> x3  0.02635229
#> 
#> $share
#> [1] 0.33241
#> 
#> $uniqueness
#> [1] 0.0050000 0.9984647 0.9993056
#> 
```
