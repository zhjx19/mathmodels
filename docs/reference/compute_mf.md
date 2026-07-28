# Build membership functions from thresholds

`compute_mf_funs` constructs a list of membership functions (one per
evaluation level) from threshold values. `compute_mf` evaluates a single
indicator value against those functions, returning a membership vector.

## Usage

``` r
compute_mf_funs(thresholds, .builder = "tri", ...)

compute_mf(x, thresholds, .builder = "tri", ...)
```

## Arguments

- thresholds:

  A numeric vector of length \\n \ge 2\\ defining the characteristic
  values (peaks / centers) of each evaluation level.

- .builder:

  A character string or function that specifies how membership functions
  are built from `thresholds`:

  `"tri"`

  :   (default) Piecewise linear (triangular + trapezoidal). First level
      decays linearly from 1 to 0 across `[th[1], th[2]]`, last level
      rises from 0 to 1 across `[th[n-1], th[n]]`, middle levels are
      isosceles triangles peaking at each threshold.

  `"gauss"`

  :   Gaussian (normal) membership. First level is a right-half Gaussian
      decaying from 1 at `th[1]`, last level is a left-half Gaussian
      rising to 1 at `th[n]`, middle levels are full Gaussians centered
      at each threshold. Pass `sigma` to control spread (default =
      `mean(diff(thresholds)) / 3`).

  `"sigmoid"`

  :   Sigmoid-based membership. First level uses a decreasing sigmoid,
      last level an increasing sigmoid, middle levels use
      difference-of-sigmoids (bell-shaped). Pass `slope` to control
      steepness (default = `5 / mean(diff(thresholds))`).

  User-supplied function

  :   A function with signature `function(thresholds, ...)` that returns
      a list of \\n\\ functions, each accepting a numeric vector `x` and
      returning membership values in \\\[0, 1\]\\.

- ...:

  Additional arguments passed to the builder (e.g. `sigma` for
  `"gauss"`, `slope` for `"sigmoid"`, or forwarded to a custom builder
  function).

- x:

  Numeric vector, input values for which to compute membership.

## Value

- `compute_mf_funs`: A list of \\n\\ functions, one per level.

- `compute_mf`: A numeric vector of length \\n\\ with membership degrees
  for each level.

## Examples

``` r
# Triangular membership (default)
th = c(0.05, 0.15, 0.25, 0.5)
compute_mf(0.07, th)
#> [1] 0.8 0.2 0.0 0.0

# Gaussian membership with custom sigma
compute_mf(0.07, th, .builder = "gauss", sigma = 0.05)
#> [1] 9.231163e-01 2.780373e-01 1.533811e-03 8.705427e-17

# Sigmoid membership with custom slope
compute_mf(0.07, th, .builder = "sigmoid", slope = 20)
#> [1] 0.4013123399 0.1883274470 0.0366957707 0.0001840719

# Custom builder: exponential decay
exp_builder = function(thresholds, rate = 1) {
  n = length(thresholds)
  lapply(seq_len(n), function(i) {
    force(i)
    function(x) exp(-rate * abs(x - thresholds[i]))
  })
}
compute_mf(0.07, th, .builder = exp_builder, rate = 10)
#> [1] 0.81873075 0.44932896 0.16529889 0.01356856

if (FALSE) { # \dontrun{
# Visualise all levels for a given builder
mfs = compute_mf_funs(th, .builder = "gauss", sigma = 0.05)
plots = lapply(mfs, \(f) plot_mf(f, xlim = c(0, 0.6)))
gridExtra::grid.arrange(grobs = plots, nrow = 2)
} # }
```
