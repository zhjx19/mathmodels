# Plot compartment trajectories

Draws one line per compartment over time. Use `compartments` to select a
subset of the state variables; the default is to plot all columns except
`time`.

## Usage

``` r
plot_compartments(df, compartments = NULL)
```

## Arguments

- df:

  A data frame returned by
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  or any built-in model function.

- compartments:

  Character vector of compartment names to plot. When `NULL` (the
  default) every column except `time` is plotted.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_compartments(sir)

plot_compartments(sir, compartments = c("S", "I"))

```
