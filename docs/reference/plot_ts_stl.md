# STL Decomposition Plot

STL Decomposition Plot

## Usage

``` r
plot_ts_stl(stl_result, title = "STL Decomposition")
```

## Arguments

- stl_result:

  Result from
  [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md).

- title:

  Plot title.

## Value

A patchwork plot.

## Examples

``` r
plot_ts_stl(ts_stl(as_ts_df(AirPassengers)))
```
