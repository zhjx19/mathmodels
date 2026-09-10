# Simple Missing Value Imputation

Fills missing values in the selected columns with a simple statistic:
group mean, median, or mode (with `.by` slice grouping), row-order
linear or cubic-spline interpolation (reusing
[`interp_linear()`](https://zhjx19.github.io/mathmodels/reference/interp_linear.md)
/
[`interp_spline()`](https://zhjx19.github.io/mathmodels/reference/interp_spline.md)),
or a constant. Returns a data frame of the same shape; the changed cells
are recorded in `attr(result, "changes")`.

## Usage

``` r
impute(
  data,
  .cols,
  method = c("mean", "median", "mode", "linear", "spline", "constant"),
  .by = NULL,
  value = NULL
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to impute.

- method:

  `"mean"`, `"median"`, `"mode"`, `"linear"`, `"spline"`, or
  `"constant"`.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; the statistic is computed within each slice.
  Supported for `"mean"`, `"median"`, and `"mode"`.

- value:

  Constant fill value (required for `method = "constant"`).

## Value

The data frame with missing values filled, with attribute `changes`
(tibble: `variable`, `row`, `old`, `new`).

## Details

`"linear"` / `"spline"` interpolate over row order, i.e. they are meant
for ordered (time-series-like) data. Missing entries at the two ends
cannot be interpolated and remain `NA`. Fully missing columns stay `NA`
with a warning. `"mode"` works for any column type; `"mean"` and
`"median"` require numeric columns.

## Examples

``` r
d = data.frame(
  g = rep(c("a", "b"), each = 3),
  x = c(1, NA, 3, 10, NA, 12)
)
impute(d, .cols = x, method = "mean", .by = g)
#>   g  x
#> 1 a  1
#> 2 a  2
#> 3 a  3
#> 4 b 10
#> 5 b 11
#> 6 b 12
impute(d, .cols = x, method = "constant", value = 0)
#>   g  x
#> 1 a  1
#> 2 a  0
#> 3 a  3
#> 4 b 10
#> 5 b  0
#> 6 b 12
```
