# Pivot Table

Builds a pivot table: aggregates `values` by the cross of `index` rows
and `columns`, then spreads the aggregates into a wide table (one column
per column level x value).

## Usage

``` r
pivot_table(
  data,
  index,
  columns,
  values,
  aggfunc = mean,
  fill_value = NA,
  names_sep = "_"
)
```

## Arguments

- data:

  A data frame.

- index:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Row identifier column(s).

- columns:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Column factor(s) to spread over.

- values:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric value column(s) to aggregate.

- aggfunc:

  Aggregation function (default `mean`).

- fill_value:

  Value used for empty cells (default `NA`).

- names_sep:

  Separator between column names and value names when several `values`
  are spread.

## Value

A tibble in wide form.

## Examples

``` r
pivot_table(ToothGrowth, index = supp, columns = dose, values = len)
#> # A tibble: 2 × 4
#>   supp  len_0.5 len_1 len_2
#>   <fct>   <dbl> <dbl> <dbl>
#> 1 VC       7.98  16.8  26.1
#> 2 OJ      13.2   22.7  26.1
pivot_table(ToothGrowth, index = supp, columns = dose, values = len,
  aggfunc = length)
#> # A tibble: 2 × 4
#>   supp  len_0.5 len_1 len_2
#>   <fct>   <int> <int> <int>
#> 1 VC         10    10    10
#> 2 OJ         10    10    10
```
