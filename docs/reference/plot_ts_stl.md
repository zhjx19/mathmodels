# Plot STL Decomposition Components

Stacked four-panel plot: observed, trend, seasonal, remainder.

## Usage

``` r
plot_ts_stl(stl_result, title = "STL Decomposition")
```

## Arguments

- stl_result:

  Result list from
  [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md).

- title:

  Character. Plot title.

## Value

A `patchwork` composite `ggplot`.

## Examples

``` r
data(AirPassengers)
res = ts_stl(AirPassengers)
plot_ts_stl(res)

```
