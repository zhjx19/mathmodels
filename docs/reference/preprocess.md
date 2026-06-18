# Preprocessing Functions for Data Normalization and Standardization

A collection of functions to preprocess numeric data, including
standardization, L2 norm normalization, Min-Max scaling, centered-type
normalization, interval-type normalization, extreme-value-based
normalization, initial-value-based normalization, mean-based
normalization, and negative-to-positive transformation. These functions
transform a numeric vector to a standardized or normalized scale,
suitable for various indicator types (positive, negative, centered,
interval-based, or extreme-based).

## Usage

``` r
standardize(x, center = TRUE, scale = TRUE)

normalize(x)

rescale(x, type = "+", a = 0, b = 1)

rescale_middle(x, m)

rescale_interval(x, a, b)

rescale_extreme(x, type = "+")

rescale_initial(x, type = "+")

rescale_mean(x)

to_positive(x, type = "minmax")
```

## Arguments

- x:

  Numeric vector to be preprocessed.

- center:

  Logical or numeric scalar, passed to
  [`base::scale`](https://rdrr.io/r/base/scale.html) for centering (for
  `standardize`). Default is `TRUE`.

- scale:

  Logical or numeric scalar, passed to
  [`base::scale`](https://rdrr.io/r/base/scale.html) for scaling (for
  `standardize`). Default is `TRUE`.

- type:

  Character scalar specifying the transformation direction or method:

  "+"

  :   Positive direction (larger values are better, for `rescale`,
      `rescale_extreme` and `rescale_initial`).

  "-"

  :   Negative direction (smaller values are better, for `rescale`
      `rescale_extreme` and `rescale_initial`).

  "minmax"

  :   Min-max transformation (for `to_positive`).

  "reciprocal"

  :   Reciprocal transformation (for `to_positive`).

- a:

  Numeric scalar, lower bound of the output range or optimal interval
  (for `rescale` and `rescale_interval`).

- b:

  Numeric scalar, upper bound of the output range or optimal interval
  (for `rescale` and `rescale_interval`).

- m:

  Numeric scalar, optimal value for centered-type normalization (for
  `rescale_middle`).

## Value

A numeric vector of the same length as `x`, transformed as follows:

- `standardize`: Standardized values (mean = 0, sd = 1).

- `normalize`: L2 norm normalized values (Euclidean norm, unit length).

- `rescale`: Min-Max scaled values in `[a, b]` (default `[0, 1]`).

- `rescale_middle`: Centered-type normalized values in `[0, 1]`, where 1
  indicates `x = m`.

- `rescale_interval`: Interval-type normalized values in `[0, 1]`, where
  1 indicates `x` in `[a, b]`.

- `rescale_extreme`: Extreme-based normalized values using `min(x)/x`
  (positive) or `x/max(x)` (negative).

- `rescale_initial`: Initial-based normalized values using `x/x[1]` or
  `x[1]/x`.

- `rescale_mean`: Mean-based normalized values using `x/mean(x)`.

- `to_positive`: Transformed values converting negative indicators to
  positive using min-max or reciprocal transformation.

## Details

These functions support various preprocessing needs in data analysis:

- `standardize`: Applies Z-score standardization (mean = 0, sd = 1),
  ideal for equalizing variances or normally distributed data.

- `normalize`: Scales the vector to unit length by dividing by its L2
  (Euclidean) norm, useful for machine learning or similarity
  calculations.

- `rescale`: Performs Min-Max scaling to a specified range (default
  `[0, 1]`), supporting positive or negative indicators.

- `rescale_middle`: Normalizes centered-type indicators, where values
  closer to an optimal value `m` are better, mapping to `[0, 1]`.

- `rescale_interval`: Normalizes interval-type indicators, where values
  within `[a, b]` are optimal, mapping to `[0, 1]`.

- `rescale_extreme`: Normalizes using extreme values: `min(x)/x` for
  positive indicators or `x/max(x)` for negative indicators, often used
  in grey relational analysis.

- `rescale_initial`: Normalizes by dividing by the first value (`x/x[1]`
  or `x[1]/x`), commonly used in grey relational analysis.

- `rescale_mean`: Normalizes by dividing by the mean (`x/mean(x)`),
  commonly used in grey relational analysis.

- `to_positive`: Converts negative indicators to positive using either
  min-max (`max(x) - x`) or reciprocal (`1/x`) transformation.

Missing values (`NA`) are preserved in the output. For `rescale_initial`
and `rescale_mean`, the initial value or mean must be non-zero,
respectively.

## Examples

``` r
# Standardization
x = c(4, 1, NA, 5, 8)
standardize(x)
#> [1] -0.1732051 -1.2124356         NA  0.1732051  1.2124356

# L2 norm normalization
normalize(x)
#> [1] 0.38851434 0.09712859         NA 0.48564293 0.77702869

# Min-Max normalization (positive direction)
rescale(x)                # Scale to \code{[0, 1]}
#> [1] 0.4285714 0.0000000        NA 0.5714286 1.0000000
rescale(x, type = "-", a = 0.002, b = 0.996)  # Reverse scaling
#> [1] 0.570 0.996    NA 0.428 0.002

# Negative-to-positive transformation
to_positive(x)                       # Min-max transformation
#> [1]  4  7 NA  3  0
to_positive(x, type = "reciprocal")  # Reciprocal transformation
#> [1] 0.250 1.000    NA 0.200 0.125

# Centered-type normalization
PH = 6:9
rescale_middle(PH, 7)
#> [1] 0.5 1.0 0.5 0.0

# Interval-type normalization
Temp = c(35.2, 35.8, 36.6, 37.1, 37.8, 38.4)
rescale_interval(Temp, 36, 37)
#> [1] 0.4285714 0.8571429 1.0000000 0.9285714 0.4285714 0.0000000

# Extreme-based normalization
rescale_extreme(x)         # min(x)/x
#> [1] 0.500 0.125    NA 0.625 1.000
rescale_extreme(x, "-")    # x/max(x)
#> [1] 0.250 1.000    NA 0.200 0.125

# Initial-based normalization
rescale_initial(x)
#> [1] 1.00 0.25   NA 1.25 2.00

# Mean-based normalization
rescale_mean(x)
#> [1] 0.8888889 0.2222222        NA 1.1111111 1.7777778
```
