# Linear and Quadratic Discriminant Analysis

Fits LDA or QDA ([`MASS::lda`](https://rdrr.io/pkg/MASS/man/lda.html) /
[`MASS::qda`](https://rdrr.io/pkg/MASS/man/qda.html)) on selected
numeric columns, reports training predictions with a confusion matrix
and accuracy, and can classify new observations.

## Usage

``` r
mv_discrim(
  data,
  group,
  cols,
  method = c("lda", "qda"),
  prior = NULL,
  newdata = NULL
)
```

## Arguments

- data:

  A data frame.

- group:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The single grouping column (coerced to a factor; must not contain
  `NA`).

- cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric predictor columns; must not include the grouping column
  (exclude it with `-group`).

- method:

  Either `"lda"` (default) or `"qda"`.

- prior:

  Optional prior probabilities for the groups, passed to
  [`MASS::lda()`](https://rdrr.io/pkg/MASS/man/lda.html) /
  [`MASS::qda()`](https://rdrr.io/pkg/MASS/man/qda.html).

- newdata:

  Optional data frame with the same predictor columns to classify; see
  `new_predictions` in the return value.

## Value

An object of class `mv_discrim`: a list with components `predictions`
(tibble with `.row`, `actual`, `predicted`, and `post_<level>` posterior
columns), `confusion` (long tibble with `actual`, `predicted`, `n`),
`scores` (tibble with `.row` and `LD1`... — LDA only, `NULL` for QDA),
`structure` (tibble with `variable` and `LD1`... for LDA, or per-group
means for QDA), `new_predictions` (tibble with `.row` — the row number
within `newdata` — plus `predicted` and posterior columns; rows of
`newdata` with missing predictor values are skipped), and `meta`.

## Examples

``` r
r = mv_discrim(iris, Species, Sepal.Length:Petal.Width)
r$meta$accuracy
#> [1] 0.98
r$confusion
#> # A tibble: 9 × 3
#>   actual     predicted      n
#>   <fct>      <fct>      <int>
#> 1 setosa     setosa        50
#> 2 versicolor setosa         0
#> 3 virginica  setosa         0
#> 4 setosa     versicolor     0
#> 5 versicolor versicolor    48
#> 6 virginica  versicolor     1
#> 7 setosa     virginica      0
#> 8 versicolor virginica      2
#> 9 virginica  virginica     49
```
