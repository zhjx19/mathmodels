# Compute epidemiological decision metrics

This function extracts key epidemiological indicators from compartmental
model outputs for decision support.

## Usage

``` r
epidemic_metrics(df, params, N = NULL, threshold = NULL)
```

## Arguments

- df:

  Data frame with time series from ODE solver. Must contain at least
  columns: time, S, I.

- params:

  Named vector or list with:

  beta

  :   Transmission rate

  gamma

  :   Recovery rate

- N:

  Population size (optional). If NULL, inferred from initial state.

- threshold:

  Healthcare or policy threshold for I(t)

## Value

A list containing:

- summary:

  Key scalar metrics

- trajectory_metrics:

  Time-varying derived indicators
