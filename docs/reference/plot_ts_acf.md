# Autocorrelation Plot

Autocorrelation Plot

## Usage

``` r
plot_ts_acf(x, max_lag = 40L, title = "ACF")
```

## Arguments

- x:

  A complete `ts_df` with no missing values.

- max_lag:

  Maximum lag.

- title:

  Plot title.

## Value

A ggplot object.

## Examples

``` r
plot_ts_acf(as_ts_df(log(AirPassengers)))
```
