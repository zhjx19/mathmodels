# Plot cumulative infections

Computes cumulative infections as \\S(0) - S(t)\\ and plots the result
over time. Requires a column named `"S"`.

## Usage

``` r
plot_cumulative_infection(df)
```

## Arguments

- df:

  A data frame with columns `time` and `S`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_cumulative_infection(sir)

```
