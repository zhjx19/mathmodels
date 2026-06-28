# Dual-Axis Plot for SARIMA-GARCH: Mean + Conditional Volatility

Upper panel shows observed values with mean fitted line; lower panel
shows the conditional standard deviation from the GARCH component.

## Usage

``` r
plot_ts_sarima_garch(
  sg_result,
  title = "SARIMA-GARCH: Mean & Volatility",
  y_lab = "Value"
)
```

## Arguments

- sg_result:

  Result list from
  [`ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima_garch.md).

- title:

  Character. Overall plot title.

- y_lab:

  Character. y-axis label for upper panel.

## Value

A `patchwork` composite `ggplot` object.

## Examples

``` r
data(AirPassengers)
res = ts_sarima_garch(log(AirPassengers))
plot_ts_sarima_garch(res)

```
