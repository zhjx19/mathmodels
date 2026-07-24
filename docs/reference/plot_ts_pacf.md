# Partial Autocorrelation Plot

Partial Autocorrelation Plot

## Usage

``` r
plot_ts_pacf(x, max_lag = 40L, title = "PACF")
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
plot_ts_pacf(as_ts_df(log(AirPassengers)))
```
