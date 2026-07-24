# Stationarity Tests for a Time Series

Runs ADF and KPSS tests on a complete `ts_df`.

## Usage

``` r
ts_test(x, adf_lags = NULL)
```

## Arguments

- x:

  A complete `ts_df` with no missing values.

- adf_lags:

  Optional lag for the ADF test.

## Value

A tibble with test results.

## Examples

``` r
ts_test(as_ts_df(log(AirPassengers)))
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo 
#> # A tibble: 2 × 5
#>   test  null_hypothesis                    statistic p_value conclusion    
#>   <chr> <chr>                                  <dbl>   <dbl> <chr>         
#> 1 ADF   Unit root present (non-stationary)     -6.42    0.01 Stationary    
#> 2 KPSS  Series is stationary                    2.83    0.01 Non-stationary
```
