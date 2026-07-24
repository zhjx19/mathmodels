# Plot a Time Series

Plot a Time Series

## Usage

``` r
plot_ts(x, title = "Time Series", x_lab = "Time", y_lab = "Value")
```

## Arguments

- x:

  A `ts_df`.

- title:

  Plot title.

- x_lab, y_lab:

  Axis labels.

## Value

A ggplot object.

## Examples

``` r
plot_ts(as_ts_df(AirPassengers))
```
