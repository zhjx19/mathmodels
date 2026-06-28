# Transform a Time Series for Stationarity

Applies variance-stabilising transformation (Box-Cox / log) and/or
differencing, returning the transformed series together with all
parameters needed to back-transform forecasts.

## Usage

``` r
ts_transform(
  x,
  method = c("none", "log", "boxcox"),
  lambda = NULL,
  diff = 0L,
  seasonal_diff = 0L,
  frequency = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  A numeric vector or `ts` object.

- method:

  Character. Variance transformation: `"none"` (default), `"log"`,
  `"boxcox"` (auto-selects lambda via
  [`forecast::BoxCox.lambda()`](https://pkg.robjhyndman.com/forecast/reference/BoxCox.lambda.html)).

- lambda:

  Numeric. Fixed Box-Cox lambda. Ignored unless `method = "boxcox"`. If
  `NULL` (default), lambda is estimated.

- diff:

  Integer. Number of regular differences (default `0`).

- seasonal_diff:

  Integer. Number of seasonal differences (default `0`).

- frequency:

  Integer. Series frequency, required when `x` is a plain numeric vector
  and `seasonal_diff > 0`.

- verbose:

  Logical. Print transformation summary (default `TRUE`).

## Value

A named list:

- transformed:

  Transformed `ts` object ready for modelling.

- params:

  Named list with `method`, `lambda`, `diff`, `seasonal_diff`,
  `frequency` — everything needed by
  [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md).

- summary:

  One-row tibble describing each step applied.

## See also

[`ts_back_transform`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md)

## Examples

``` r
data(AirPassengers)
res = ts_transform(AirPassengers, method = "log", diff = 1,
                    seasonal_diff = 1)
#> 
#> === ts_transform() ===
#> Steps: log(x) -> seasonal_diff(lag=12) -> diff(lag=1) 
#> Length: 144 -> 131 
res$summary
#> # A tibble: 1 × 4
#>   original_length transformed_length steps_applied                        lambda
#>             <int>              <int> <chr>                                 <dbl>
#> 1             144                131 log(x) -> seasonal_diff(lag=12) -> …     NA
plot_ts(res$transformed, title = "Transformed AirPassengers")

```
