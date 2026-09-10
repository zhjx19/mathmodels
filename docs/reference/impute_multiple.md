# Multiple Imputation via FCS

Generates `m` completed data sets by fully conditional specification
(chained equations). Numeric columns use proper Bayesian linear-model
draws with predicted-mean matching (`"pmm"`, default) or normal noise
(`"norm"`); binary and polytomous columns use parameter-perturbed
logistic (`"logreg"`) / multinomial (`"polyreg"`) models — approximately
proper, following the same spirit as `mice`.

## Usage

``` r
impute_multiple(
  data,
  .cols = dplyr::everything(),
  m = 5,
  maxit = 5,
  seed = NULL,
  methods = NULL
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to impute; defaults to all columns. Columns without missing
  values are left as is.

- m:

  Number of imputations (at least 2).

- maxit:

  Number of chain iterations per imputation.

- seed:

  Random seed for reproducibility.

- methods:

  Optional named vector/list overriding the automatic method choice per
  column, e.g. `methods = c(x = "norm")`.

## Value

An object of class `mi_impute`: a list of `m` completed data frames.
Attributes: `methods`, `iterations`, `trace` (tibble with imputation /
iteration / variable / mean / sd of the imputed cells; mean/sd are NA
for categorical columns), `seed`, `n_imputed`.

## Details

Logical and character columns are treated as categorical and come back
as factors. Convergence can be inspected via the `trace` attribute
(plotting is planned for a later version).

## Examples

``` r
set.seed(1)
d = data.frame(x = rnorm(30), z = rnorm(30),
               g = factor(sample(c("a", "b"), 30, replace = TRUE)))
d$x[c(2, 8)] = NA
d$g[5] = NA
imps = impute_multiple(d, m = 3, maxit = 3, seed = 11)
imps
#> Multiple imputation (mi_impute): 3 imputation(s) of 2 column(s)
#>   methods: x=pmm, g=logreg 
#>   cells imputed: x=2, g=1 
#>   iterations: 3 
```
