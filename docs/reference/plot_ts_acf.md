# ACF and PACF Plots

ACF and PACF Plots

## Usage

``` r
plot_ts_acf(
  x,
  max_lag = 40L,
  type = c("acf", "pacf", "both"),
  title = NULL,
  diff = 0L
)
```

## Arguments

- x:

  A numeric vector, `ts` object, or tidy tibble/data.frame with a
  `value` column.

- max_lag:

  Integer. Maximum lag (default 40).

- type:

  Character. `"acf"` (default), `"pacf"`, or `"both"` (returns a
  patchwork panel).

- title:

  Character. Plot title.

- diff:

  Integer. Apply [`diff()`](https://rdrr.io/r/base/diff.html) this many
  times before plotting (default 0).

## Value

A `ggplot` or `patchwork` object.

## Examples

``` r
data(AirPassengers)
plot_ts_acf(log(AirPassengers), type = "both")

```
