# Model-Based Missing Value Imputation

Unified imputation framework on top of three models: k-nearest-neighbors
(`knn`, hand-written), linear regression (`lm`,
[`stats::lm`](https://rdrr.io/r/stats/lm.html)), and regression trees
(`tree`, [`rpart::rpart`](https://rdrr.io/pkg/rpart/man/rpart.html)).
For every numeric column with missing values, a model is trained on the
complete rows (target column and all other numeric columns as features)
and the missing cells are replaced by the predictions. Returns a data
frame of the same shape with the changed cells recorded in
`attr(result, "changes")` and per-column details in
`attr(result, "details")`.

## Usage

``` r
impute_model(
  data,
  .cols,
  model = c("knn", "lm", "tree"),
  .by = NULL,
  k = 5,
  ...
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric target columns to impute. Features are all other numeric
  columns of `data`.

- model:

  `"knn"`, `"lm"`, or `"tree"`.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; models are trained separately within each
  slice.

- k:

  Number of neighbors for `model = "knn"`.

- ...:

  Additional arguments passed to
  [`rpart::rpart()`](https://rdrr.io/pkg/rpart/man/rpart.html) for
  `model = "tree"` (e.g. `control = rpart::rpart.control(cp = 0.01)`).

## Value

The data frame with missing values filled, with attributes `changes`
(tibble: `variable`, `row`, `old`, `new`) and `details` (tibble:
`variable`, `model`, `n_train`, `n_imputed`).

## Details

Rows whose features still contain missing values cannot be scored by
`"lm"` / `"tree"` predictions and keep their `NA` (a joint iterative
scheme is planned for a later version). `"knn"` handles rows with
partially missing features by computing distances on the features
observed for the target row; rows with all features missing keep `NA`.

## Examples

``` r
set.seed(1)
d = data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
d$x[c(3, 10)] = NA
impute_model(d, .cols = x, model = "knn", k = 5)
#>              x           y            z
#> 1  -0.62645381  1.35867955  2.401617761
#> 2   0.18364332 -0.10278773 -0.039240003
#> 3   0.18097942  0.38767161  0.689739362
#> 4   1.59528080 -0.05380504  0.028002159
#> 5   0.32950777 -1.37705956 -0.743273209
#> 6  -0.82046838 -0.41499456  0.188792300
#> 7   0.48742905 -0.39428995 -1.804958629
#> 8   0.73832471 -0.05931340  1.465554862
#> 9   0.57578135  1.10002537  0.153253338
#> 10 -0.11812004  0.76317575  2.172611670
#> 11  1.51178117 -0.16452360  0.475509529
#> 12  0.38984324 -0.25336168 -0.709946431
#> 13 -0.62124058  0.69696338  0.610726353
#> 14 -2.21469989  0.55666320 -0.934097632
#> 15  1.12493092 -0.68875569 -1.253633400
#> 16 -0.04493361 -0.70749516  0.291446236
#> 17 -0.01619026  0.36458196 -0.443291873
#> 18  0.94383621  0.76853292  0.001105352
#> 19  0.82122120 -0.11234621  0.074341324
#> 20  0.59390132  0.88110773 -0.589520946
#> 21  0.91897737  0.39810588 -0.568668733
#> 22  0.78213630 -0.61202639 -0.135178615
#> 23  0.07456498  0.34111969  1.178086997
#> 24 -1.98935170 -1.12936310 -1.523566800
#> 25  0.61982575  1.43302370  0.593946188
#> 26 -0.05612874  1.98039990  0.332950371
#> 27 -0.15579551 -0.36722148  1.063099837
#> 28 -1.47075238 -1.04413463 -0.304183924
#> 29 -0.47815006  0.56971963  0.370018810
#> 30  0.41794156 -0.13505460  0.267098791
impute_model(d, .cols = x, model = "lm")
#>              x           y            z
#> 1  -0.62645381  1.35867955  2.401617761
#> 2   0.18364332 -0.10278773 -0.039240003
#> 3   0.18738733  0.38767161  0.689739362
#> 4   1.59528080 -0.05380504  0.028002159
#> 5   0.32950777 -1.37705956 -0.743273209
#> 6  -0.82046838 -0.41499456  0.188792300
#> 7   0.48742905 -0.39428995 -1.804958629
#> 8   0.73832471 -0.05931340  1.465554862
#> 9   0.57578135  1.10002537  0.153253338
#> 10  0.30097327  0.76317575  2.172611670
#> 11  1.51178117 -0.16452360  0.475509529
#> 12  0.38984324 -0.25336168 -0.709946431
#> 13 -0.62124058  0.69696338  0.610726353
#> 14 -2.21469989  0.55666320 -0.934097632
#> 15  1.12493092 -0.68875569 -1.253633400
#> 16 -0.04493361 -0.70749516  0.291446236
#> 17 -0.01619026  0.36458196 -0.443291873
#> 18  0.94383621  0.76853292  0.001105352
#> 19  0.82122120 -0.11234621  0.074341324
#> 20  0.59390132  0.88110773 -0.589520946
#> 21  0.91897737  0.39810588 -0.568668733
#> 22  0.78213630 -0.61202639 -0.135178615
#> 23  0.07456498  0.34111969  1.178086997
#> 24 -1.98935170 -1.12936310 -1.523566800
#> 25  0.61982575  1.43302370  0.593946188
#> 26 -0.05612874  1.98039990  0.332950371
#> 27 -0.15579551 -0.36722148  1.063099837
#> 28 -1.47075238 -1.04413463 -0.304183924
#> 29 -0.47815006  0.56971963  0.370018810
#> 30  0.41794156 -0.13505460  0.267098791
```
