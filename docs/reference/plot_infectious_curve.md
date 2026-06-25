# Plot infectious population

Draws a single line showing the number of infectious individuals over
time. Requires a column named `"I"`.

## Usage

``` r
plot_infectious_curve(df)
```

## Arguments

- df:

  A data frame with columns `time` and `I`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_infectious_curve(sir)

```
