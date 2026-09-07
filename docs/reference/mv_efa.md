# Exploratory Factor Analysis

Performs an exploratory factor analysis on selected numeric columns via
maximum likelihood
([`stats::factanal`](https://rdrr.io/r/stats/factanal.html)) on the
standardized data, and returns tidy loadings, uniquenesses, variance
shares, and (optionally) factor scores.

## Usage

``` r
mv_efa(
  data,
  cols,
  factors,
  rotation = "varimax",
  scores = c("regression", "none")
)
```

## Arguments

- data:

  A data frame.

- cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns to analyse.

- factors:

  Number of factors to extract.

- rotation:

  Character. Rotation passed to `factanal`: one of `"varimax"`
  (default), `"promax"`, `"none"`, or the name of a rotation function.

- scores:

  Either `"regression"` (default) to compute factor scores, or `"none"`.

## Value

An object of class `mv_efa`: a list with components `loadings` (tibble
with `variable`, `Factor1`..., and `uniqueness`), `uniquenesses`,
`variance` (tibble with `factor`, `ss_loadings`, `prop_var`, `cum_var`),
`scores` (tibble with `.row` and `Factor1`..., or `NULL`), and `meta`.

## Examples

``` r
r = mv_efa(mtcars, cyl:carb, factors = 2)
r$loadings
#> # A tibble: 10 × 4
#>    variable Factor1 Factor2 uniqueness
#>    <chr>      <dbl>   <dbl>      <dbl>
#>  1 cyl      -0.633   0.740      0.0522
#>  2 disp     -0.728   0.610      0.0974
#>  3 hp       -0.330   0.854      0.162 
#>  4 drat      0.810  -0.231      0.291 
#>  5 wt       -0.795   0.411      0.198 
#>  6 qsec     -0.163  -0.915      0.137 
#>  7 vs        0.292  -0.821      0.240 
#>  8 am        0.903   0.0783     0.179 
#>  9 gear      0.869   0.114      0.232 
#> 10 carb      0.0493  0.768      0.407 
```
