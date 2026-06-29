# Extract key epidemic metrics from ODE simulation output

This function computes core epidemiological summary statistics from
deterministic compartmental epidemic models (e.g., SIR, SIS, SEIR).

## Usage

``` r
epi_metrics(data, beta, gamma, N = NULL)
```

## Arguments

- data:

  A data frame returned by `model_*()` functions. Must contain a column
  `time`, and at least `S` and `I`. Additional compartments such as `E`
  or `R` are allowed.

- beta:

  Transmission rate (numeric).

- gamma:

  Recovery / removal rate (numeric).

- N:

  Total population (numeric). If `NULL` (default), it is computed as the
  sum of all compartment values at the first time point.

## Value

A named list with:

- R0:

  Basic reproduction number.

- peak_infection:

  Maximum number of infectious individuals.

- peak_time:

  Time at which infectious peak occurs.

- attack_rate:

  Proportion of population that transitioned out of S; defined only for
  models containing an R compartment (e.g., SIR/SEIR).

## Details

The input data frame is expected to be the output of `model_*()`
functions, containing a time column and one or more compartment columns
(e.g., S, I, R, E). All compartment values are assumed to be in absolute
population counts, and the total population is assumed to be conserved
unless otherwise specified.

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
