# Little's MCAR Test

Performs Little's (1988) test of Missing Completely At Random (the
implementation follows
[`naniar::mcar_test()`](https://naniar.njtierney.com/reference/mcar_test.html)
/ Hawkins' refinement): the mean vector and covariance matrix are first
estimated under the MCAR assumption by maximum likelihood (EM
algorithm), then each missingness pattern's mean (on its observed
columns) is compared with the overall mean via a Mahalanobis distance
weighted by the pattern size. The statistic is referred to a chi-square
distribution with `df = sum(observed columns per pattern) - p`.

## Usage

``` r
na_mcar_test(data, .cols = dplyr::everything(), maxit = 200, tol = 1e-08)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns to test; defaults to all columns.

- maxit:

  Maximum EM iterations.

- tol:

  EM convergence tolerance.

## Value

A one-row tibble: `chi2`, `df`, `p_value`, `n_patterns`, `n_used`,
`p_missing_cells`. Attributes: `patterns` — a tibble with one row per
missingness pattern (missing columns, number of rows).

## Details

The EM estimate requires an invertible covariance matrix; fully missing
columns are dropped with a warning before fitting. With no missing
values (or a degenerate pattern structure) the test is not applicable
and an error is raised.

## Examples

``` r
set.seed(1)
d = data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
d$y[1:10] = NA
na_mcar_test(d)
#> # A tibble: 1 × 6
#>    chi2    df p_value n_patterns n_used p_missing_cells
#> * <dbl> <int>   <dbl>      <int>  <int>           <dbl>
#> 1  3.95     2   0.139          2     50          0.0667
```
