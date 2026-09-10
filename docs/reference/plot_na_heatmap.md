# Missingness Heatmap

Draws a row-by-column missingness heatmap: each cell shows whether the
entry is observed (grey) or missing (red). Rows keep their original
order so the plot stays aligned with the data.

## Usage

``` r
plot_na_heatmap(data, .cols)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to display.

## Value

A `ggplot` object.

## Examples

``` r
d = data.frame(x = c(1, NA, 3, 4), y = c(NA, 2, NA, 5))
plot_na_heatmap(d)
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.
```
