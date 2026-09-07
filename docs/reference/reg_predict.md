# Predictions

Generates point predictions and confidence/ prediction intervals for new
data.

## Usage

``` r
reg_predict(model_result, n_new = 10, level = 0.95, newdata = NULL)
```

## Arguments

- model_result:

  A fitted result from any `reg_*()` function.

- n_new:

  Number of new observations to generate when `newdata` is not supplied.
  Predictor values are drawn from the empirical distributions of the
  input columns (normal approximation by mean/sd for numeric columns;
  levels sampled with observed proportions for factor columns).

- level:

  Confidence/interval level. Defaults to 0.95.

- newdata:

  Optional data frame of new predictor values. When supplied,
  predictions are made for exactly these rows and `n_new` is ignored.

## Value

A named list with:

- predictions:

  Tibble with predicted values.

- confidence_interval:

  Tibble with lower and upper bounds of the confidence interval.

- prediction_interval:

  Tibble with lower and upper bounds of the prediction interval (linear
  models only).

- newdata:

  The data frame actually used for prediction.
