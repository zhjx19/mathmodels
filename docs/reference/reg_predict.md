# Predictions

Generates point predictions and confidence/ prediction intervals for new
data.

## Usage

``` r
reg_predict(model_result, n_new = 10, level = 0.95)
```

## Arguments

- model_result:

  A fitted result from any `reg_*()` function.

- n_new:

  Number of new observations to predict. Random values drawn from the
  fitted model's residuals are used for predictor data.

- level:

  Confidence/interval level. Defaults to 0.95.

## Value

A named list with:

- predictions:

  Tibble with predicted values.

- confidence_interval:

  Tibble with lower and upper bounds.
