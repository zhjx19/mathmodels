# STL (Seasonal + Trend + Loess) Decomposition

Wraps [`stats::stl()`](https://rdrr.io/r/stats/stl.html) and returns
both a tidy long-format tibble and the raw `stl` object for further use.

## Usage

``` r
ts_stl(x, frequency = NULL, s_window = "periodic", robust = TRUE, ...)
```

## Arguments

- x:

  A `ts` object with `frequency > 1`, or a numeric vector with
  `frequency` specified.

- frequency:

  Integer. Required if `x` is a plain numeric vector.

- s_window:

  Seasonal smoothing window passed to
  [`stl()`](https://rdrr.io/r/stats/stl.html). Use `"periodic"`
  (default) for purely periodic seasonality.

- robust:

  Logical. Use robust fitting in STL (default `TRUE`).

- ...:

  Additional arguments forwarded to
  [`stats::stl()`](https://rdrr.io/r/stats/stl.html).

## Value

A named list:

- components:

  Tidy `tibble` with columns
  `index, observed, trend, seasonal, remainder`.

- strength:

  Tibble: `seasonal_strength` and `trend_strength` (Wang et al. 2006
  metrics).

- model:

  Raw `stl` object.

## Examples

``` r
data(AirPassengers)
res = ts_stl(AirPassengers)
res$components
#> # A tibble: 144 × 5
#>    index observed trend seasonal remainder
#>    <int>    <dbl> <dbl>    <dbl>     <dbl>
#>  1     1      112  123.  -16.5       5.33 
#>  2     2      118  123.  -27.3      21.9  
#>  3     3      132  124.    9.01     -0.667
#>  4     4      129  124.    0.992     3.87 
#>  5     5      121  125.    3.08     -6.70 
#>  6     6      135  125.   20.0     -10.2  
#>  7     7      148  126.   33.4     -11.2  
#>  8     8      148  126.   40.1     -18.4  
#>  9     9      136  127.   22.2     -13.1  
#> 10    10      119  128.  -13.6       4.82 
#> # ℹ 134 more rows
```
