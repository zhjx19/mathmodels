# Simple Correspondence Analysis

Builds a contingency table from two columns of a data frame and runs
simple correspondence analysis
([`MASS::corresp`](https://rdrr.io/pkg/MASS/man/corresp.html)),
returning row and column coordinates and the inertia decomposition.

## Usage

``` r
mv_corresp(data, row, col, nf = 2)
```

## Arguments

- data:

  A data frame.

- row:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The single row-factor column.

- col:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The single column-factor column.

- nf:

  Number of dimensions to compute; must be smaller than both the number
  of rows and the number of columns of the table.

## Value

An object of class `mv_corresp`: a list with components `row_coord`
(tibble with `row` and `Dim1`...), `col_coord` (tibble with `col` and
`Dim1`...), `inertia` (tibble with `dimension`, `inertia`, `prop`,
`cum`), and `meta`.

## Examples

``` r
set.seed(42)
df = data.frame(
  breed = rep(c("A", "B", "C"), each = 30),
  grade = sample(c("low", "mid", "high"), 90, TRUE)
)
r = mv_corresp(df, breed, grade)
r$inertia
#> # A tibble: 2 × 4
#>   dimension  inertia     prop   cum
#>       <int>    <dbl>    <dbl> <dbl>
#> 1         1 3.97e- 2 1   e+ 0     1
#> 2         2 3.41e-32 8.59e-31     1
```
