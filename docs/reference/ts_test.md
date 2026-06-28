# Stationarity Tests for a Time Series

Runs ADF, KPSS, and PP tests and returns a tidy summary data frame.

## Usage

``` r
ts_test(
  x,
  adf_lags = NULL,
  kpss_type = c("Level", "Trend"),
  pp_type = c("Z(t_alpha)", "Z(alpha)"),
  verbose = TRUE
)
```

## Arguments

- x:

  A numeric vector, `ts` object, or single-column data frame.

- adf_lags:

  Integer. Max lag for ADF (default: `NULL` = auto via
  `trunc((length(x)-1)^(1/3))`).

- kpss_type:

  Character. `"Level"` (default) or `"Trend"` for KPSS null hypothesis.

- pp_type:

  Character. `"Z(t_alpha)"` (default) or `"Z(alpha)"`.

- verbose:

  Logical. Print test results to console (default `TRUE`).

## Value

A `tibble` with columns:

- test:

  Test name (ADF / KPSS / PP)

- null_hypothesis:

  Plain-English description of H0

- statistic:

  Test statistic

- p_value:

  p-value (or NA when only bounds are available)

- conclusion:

  Character: "Stationary" or "Non-stationary"

## Examples

``` r
data(AirPassengers)
ts_test(AirPassengers)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo 
#> Warning: p-value smaller than printed p-value
#> Warning: p-value smaller than printed p-value
#> Warning: p-value smaller than printed p-value
#> 
#> === Stationarity Tests ===
#> # A tibble: 3 × 5
#>   test  null_hypothesis                    statistic p_value conclusion    
#>   <chr> <chr>                                  <dbl>   <dbl> <chr>         
#> 1 ADF   Unit root present (non-stationary)     -7.32    0.01 Stationary    
#> 2 KPSS  Series is stationary                    2.74    0.01 Non-stationary
#> 3 PP    Unit root present (non-stationary)     -5.05    0.01 Stationary    
#> 
#> Note: All tests use alpha = 0.05.
```
