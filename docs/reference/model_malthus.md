# Malthusian (Exponential) Growth Model

Solves the Malthusian growth equation: \$\$\frac{dN}{dt} = r \\ N\$\$

## Usage

``` r
model_malthus(init, params, times, ...)
```

## Arguments

- init:

  Named numeric vector, e.g. `c(N = 100)`.

- params:

  Named numeric vector, e.g. `c(r = 0.3)`.

  r

  :   Intrinsic growth rate (per time unit).

- times:

  Numeric vector of output times.

- ...:

  Additional arguments passed to `ode_solver`.

## Value

A `data.frame` with columns `time` and `N`.

## Details

The analytical solution is \\N(t) = N_0 \\ e^{r t}\\, which can be used
to verify numerical accuracy.

## Examples

``` r
res = model_malthus(
  init   = c(N = 100),
  params = c(r = 0.3),
  times  = seq(0, 10, by = 0.1)
)
# Compare with analytical solution
res$N_exact = 100 * exp(0.3 * res$time)
head(res)
#>   time        N  N_exact
#> 1  0.0 100.0000 100.0000
#> 2  0.1 103.0455 103.0455
#> 3  0.2 106.1838 106.1837
#> 4  0.3 109.4175 109.4174
#> 5  0.4 112.7498 112.7497
#> 6  0.5 116.1835 116.1834
```
