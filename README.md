
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mathmodels

<!-- badges: start -->
<!-- badges: end -->

The goal of mathmodels is to create a mathematical modeling toolbox with
R, implement common mathematical modeling algorithms: evaluation,
optimization, prediction, dynamics, graph theory, statistics,
intelligence, etc.

**Attention:** The package is currently only a very very early version!

## Installation

You can install the latest development version of `mathmodels` directly
from [GitHub](https://github.com/zhjx19/mathmodes) use:

``` r
library(devtools) 
devtools::install_github("zhjx19/mathmodels")
```

Or download to current path, unzip and use:

``` r
install.packages("mathmodels-master", repos=NULL, type="source")
```

## Example

This is a basic example which shows you how to solve a AHP problem:

``` r
library(mathmodels)
# pairwise comparison matrix
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
# return: Weight vector, Consistency ratio, Maximum eigenvalue, Consistency index
AHP(A)
#> $W
#> [1] 0.26360349 0.47583538 0.05381460 0.09806829 0.10867824
#> 
#> $CR
#> [1] 0.01609027
#> 
#> $Lmax
#> [1] 5.072084
#> 
#> $CI
#> [1] 0.0180211
```
