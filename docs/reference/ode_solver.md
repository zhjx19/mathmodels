# General ODE Solver

Solves a system of ordinary differential equations (ODEs) using
string-formula equations. Equation strings are parsed once at call time
and reused across all integration steps, so the solver is substantially
faster than approaches that call
[`parse()`](https://rdrr.io/r/base/parse.html) inside the derivative
function.

## Usage

``` r
ode_solver(init, times, equations, params = NULL, method = "lsoda", ...)
```

## Arguments

- init:

  A named numeric vector of initial values for all state variables.

- times:

  A numeric vector of output times.

- equations:

  A named character vector; names are variable names, values are
  derivative expressions (e.g. `c(S = "-beta*S*I")`).

- params:

  A named numeric vector or list of parameters referenced inside
  `equations`.

- method:

  Integration method passed to
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html). Default
  `"lsoda"`.

- ...:

  Additional arguments forwarded to
  [`deSolve::ode`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

A `data.frame` with a `time` column followed by one column per state
variable.

## Details

All equations are parsed with
[`parse()`](https://rdrr.io/r/base/parse.html) exactly once before the
integrator starts. A single pre-allocated environment is reused on every
derivative evaluation, avoiding repeated memory allocation and garbage
collection.

## Examples

``` r
# Built-in wrapper (one-liner)
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
head(sir)
#>   time        S        I         R
#> 1  0.0 990.0000 10.00000 0.0000000
#> 2  0.1 987.8242 12.06579 0.1100099
#> 3  0.2 985.2059 14.55138 0.2427141
#> 4  0.3 982.0584 17.53891 0.4027092
#> 5  0.4 978.2793 21.12521 0.5954865
#> 6  0.5 973.7487 25.42375 0.8275867

# Custom SEIR model with demography
seir_demo = ode_solver(
  init  = c(S = 1000, E = 1, I = 0, R = 0),
  times = seq(0, 200, by = 1),
  equations = c(
    S = "mu*N - beta*S*I - mu*S",
    E = "beta*S*I - alpha*E - mu*E",
    I = "alpha*E - gamma*I - mu*I",
    R = "gamma*I - mu*R"
  ),
  params = c(beta = 0.3, alpha = 0.2, gamma = 0.1, mu = 0.01, N = 1000)
)
tail(seir_demo)
#>     time     S        E        I        R
#> 196  195 0.385 47.60071 86.54675 865.6098
#> 197  196 0.385 47.60071 86.54675 865.6084
#> 198  197 0.385 47.60071 86.54675 865.6070
#> 199  198 0.385 47.60071 86.54675 865.6056
#> 200  199 0.385 47.60071 86.54675 865.6042
#> 201  200 0.385 47.60071 86.54675 865.6029
```
