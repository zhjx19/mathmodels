# Plot incidence curves

Calls
[`compute_incidence()`](https://zhjx19.github.io/mathmodels/reference/compute_incidence.md)
and draws one line for each incidence metric (new infections from S, new
infections from I, new deaths).

## Usage

``` r
plot_incidence(df)
```

## Arguments

- df:

  A data frame returned by
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  or any built-in model function.

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
