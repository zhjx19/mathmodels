# Plot effective reproduction number R_t

Estimates the time-varying effective reproduction number \\R_t\\
directly from compartment data. Two methods are available:

- `"mechanistic"`:

  Uses \\R_t = \beta S(t) / \gamma\\. This is the standard formula for
  an SIR-type model and is recommended when the model dynamics match the
  SIR structure.

- `"normalized"`:

  Uses \\R_t = R_0 \\ S(t) / N\\ with \\R_0 = \beta N / \gamma\\.
  Suitable when different normalisations are desired.

## Usage

``` r
plot_Rt_estimate(df, params, N = NULL, method = c("mechanistic", "normalized"))
```

## Arguments

- df:

  A data frame with columns `time` and `S`.

- params:

  A named numeric vector or list containing at least `beta` and `gamma`.

- N:

  Total population size. If `NULL` (the default), it is estimated from
  the sum of all state columns at the first time step.

- method:

  Estimation method: `"mechanistic"` (default) or `"normalized"`.

## Value

A [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html) object
with a dashed horizontal line at \\R_t = 1\\.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
plot_Rt_estimate(sir, params = c(beta = 0.002, gamma = 0.1))

```
