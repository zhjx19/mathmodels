# Linear weighted synthesis

Linear weighted synthesis

## Usage

``` r
linear_sum(data, w)
```

## Arguments

- data:

  Indicator matrix or data frame (rows = samples, columns = indicators)

- w:

  Weight vector (length must match ncol(data))

## Value

Vector of comprehensive scores

## Examples

``` r
data = data.frame(
  GDP = c(0.85, 0.72, 0.91, NA),
  Employment = c(0.78, 0.85, 0.67, 0.73),
  Environment = c(0.65, 0.72, NA, 0.81)
)
w = c(0.5, 0.3, 0.2)
linear_sum(data, w)
#> [1] 0.789 0.759    NA    NA
```
