# Prediction Plot

Returns a fitted-vs-observed line plot, with optional confidence bands.

## Usage

``` r
reg_plot_predict(model_result, type = c("fitted", "predicted"), show_ci = TRUE)
```

## Arguments

- model_result:

  A fitted result from any `reg_*()` function.

- type:

  Character. Currently supports `"fitted"` only.

- show_ci:

  Logical. Whether to overlay confidence intervals. Default `TRUE`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.
