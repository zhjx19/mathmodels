# Plot daily new infections (dI)

Computes the daily change in the infectious compartment \\\Delta I\\
directly via [`diff()`](https://rdrr.io/r/base/diff.html) and plots it
as a line chart. The peak is highlighted with a point marker.

## Usage

``` r
plot_incidence(df)
```

## Arguments

- df:

  A data frame returned by
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  or any built-in model function. Must contain columns `time` and `I`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_incidence(sir)

```
