# Missingness Bar Plot

Draws a bar plot of the missing-value percentage per column.

## Usage

``` r
plot_na_bar(data, .cols)
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
plot_na_bar(d)
```
