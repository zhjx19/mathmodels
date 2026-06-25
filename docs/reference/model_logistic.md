# Logistic Population Growth Model

Solves the logistic growth equation: \$\$\frac{dN}{dt} = r \left(1 -
\frac{N}{K}\right) N\$\$

## Usage

``` r
model_logistic(init, params, times, ...)
```

## Arguments

- init:

  Named numeric vector, e.g. `c(N = 10)`.

- params:

  Named numeric vector, e.g. `c(r = 0.5, K = 100)`.

  r

  :   Intrinsic growth rate.

  K

  :   Carrying capacity.

- times:

  Numeric vector of output times.

- ...:

  Additional arguments passed to `ode_solver`.

## Value

A `data.frame` with columns `time` and `N`.

## Details

The analytical solution is \\N(t) = K \\/\\ \bigl(1 + (K/N_0 - 1)\\e^{-r
t}\bigr)\\, which can be used to verify numerical accuracy.

## Examples

``` r
res = model_logistic(
  init   = c(N = 10),
  params = c(r = 0.5, K = 100),
  times  = seq(0, 20, by = 0.1)
)
head(res)
#>   time        N
#> 1  0.0 10.00000
#> 2  0.1 10.45909
#> 3  0.2 10.93669
#> 4  0.3 11.43331
#> 5  0.4 11.94947
#> 6  0.5 12.48563
```
