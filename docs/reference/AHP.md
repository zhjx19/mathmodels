# AHP: Analytic Hierarchy Process

AHP is a multi-criteria decision analysis method developed by Saaty,
which can also be used to determine indicator weights.

## Usage

``` r
AHP(A)
```

## Arguments

- A:

  a numeric matrix, i.e. pairwise comparison matrix

## Value

a list object that contains: w (Weight vector), CR (Consistency ratio),
Lmax (Maximum eigenvalue), CI (Consistency index)

## Examples

``` r
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
AHP(A)
#> $w
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
#> 
```
