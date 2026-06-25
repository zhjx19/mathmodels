# Compute incidence from ODE output

Estimates daily new infections and new deaths from the change in
compartment sizes between consecutive time steps.

## Usage

``` r
compute_incidence(
  df,
  infection_col = "I",
  death_col = "D",
  population_col = "S"
)
```

## Arguments

- df:

  A data frame with a `time` column and compartment columns.

- infection_col:

  Name of the infectious compartment column (default `"I"`).

- death_col:

  Name of the deceased compartment column (default `"D"`). If the column
  does not exist, death incidence is silently skipped.

- population_col:

  Name of the susceptible compartment column (default `"S"`). The
  negative difference of `S` is used as a second estimate of new
  infections.

## Value

A data frame with columns `time`, `new_infection` (from S depletion),
`new_infection_I` (from I change), and `new_death` (if `death_col`
exists). The first row of each difference column is `NA`.

## Examples

``` r
sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.1)
)
inc = compute_incidence(sir)
head(inc)
#>   time new_infection new_infection_I
#> 1  0.0            NA              NA
#> 2  0.1      2.175802        2.065792
#> 3  0.2      2.618291        2.485587
#> 4  0.3      3.147523        2.987528
#> 5  0.4      3.779080        3.586302
#> 6  0.5      4.530643        4.298543
```
