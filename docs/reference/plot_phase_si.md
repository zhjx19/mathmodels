# Phase plot S vs I

Draws a trajectory in the \\(S, I)\\ plane. Useful for visualising the
epidemic orbit and threshold behaviour.

## Usage

``` r
plot_phase_si(df)
```

## Arguments

- df:

  A data frame with columns `time`, `S`, and `I`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_phase_si(sir)

```
