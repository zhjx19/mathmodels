# PACF Plot (convenience alias)

Shortcut for `plot_ts_acf(x, type = "pacf")`.

## Usage

``` r
plot_ts_pacf(x, max_lag = 40L, title = NULL, diff = 0L)
```

## Arguments

- x:

  A numeric vector, `ts` object, or tidy tibble/data.frame with a
  `value` column.

- max_lag:

  Integer. Maximum lag (default 40).

- title:

  Character. Plot title.

- diff:

  Integer. Apply [`diff()`](https://rdrr.io/r/base/diff.html) this many
  times before plotting (default 0).
