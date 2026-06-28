# Plot a Time Series

Renders a clean line chart for one or more time series.

## Usage

``` r
plot_ts(
  x,
  title = "Time Series",
  subtitle = NULL,
  x_lab = "Time",
  y_lab = "Value",
  colour = .PALETTE$observed,
  add_points = FALSE
)
```

## Arguments

- x:

  A numeric vector, `ts` object, tidy tibble with columns
  `(index, value)`, or a named list of such objects for multi-series.

- title:

  Character. Plot title.

- subtitle:

  Character. Plot subtitle.

- x_lab:

  Character. x-axis label (default `"Time"`).

- y_lab:

  Character. y-axis label (default `"Value"`).

- colour:

  Character. Line colour (ignored for multi-series).

- add_points:

  Logical. Add point markers (default `FALSE`).

## Value

A `ggplot` object.

## Examples

``` r
data(AirPassengers)
plot_ts(AirPassengers, title = "Air Passengers", y_lab = "Thousands")

```
