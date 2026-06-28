# GARCH Volatility Plot

Four-panel overview: (1) observed returns with \\\pm 2\sigma\\ bands,
(2) conditional variance, (3) standardised residuals, (4) histogram of
standardised residuals.

## Usage

``` r
plot_ts_garch(garch_result, title = "GARCH Volatility Analysis")
```

## Arguments

- garch_result:

  Result list from
  [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md).

- title:

  Character. Plot title.

## Value

A `patchwork` composite `ggplot`.

## Examples

``` r
r = diff(log(AirPassengers))
res = ts_garch(r)
plot_ts_garch(res)

```
