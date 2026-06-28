# Extract key epidemic metrics

Returns a named list of 4 core scalar metrics derived from ODE output.

## Usage

``` r
epi_metrics(df, beta, gamma, N = NULL)
```

## Arguments

- df:

  A data frame with columns `time`, `S`, `I`.

- beta:

  Transmission rate (numeric).

- gamma:

  Recovery / removal rate (numeric).

- N:

  Total population (numeric). If `NULL` (default), computed as the sum
  of compartment values at the first time point. `attack_rate`.

## Value

A named list with components: `R0` (basic reproduction number),
`peak_infection` (maximum number of infectious individuals), `peak_time`
(time at which peak occurs), `attack_rate` (proportion of susceptible
that became infected).

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
epi_metrics(sir, beta = 0.002, gamma = 0.1)
#> $R0
#> [1] 20
#> 
#> $peak_infection
#> [1] 800.5652
#> 
#> $peak_time
#> [1] 4.1
#> 
#> $attack_rate
#> [1] 0.99
#> 
```
