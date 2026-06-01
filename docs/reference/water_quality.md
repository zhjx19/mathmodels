# Water Quality Dataset

A dataset containing water quality evaluation metrics for 20 rivers,
including dissolved oxygen (O2, positive indicator), pH value (PH,
centered indicator), total bacteria count (germ, negative indicator),
and plant nutrient content (nutrient, interval indicator with optimal
range 10-20). This dataset is suitable for multi-criteria decision
analysis, such as weight calculation and fuzzy comprehensive evaluation
in the `mathmodels` package.

## Usage

``` r
water_quality
```

## Format

A data frame with 20 rows and 5 columns:

- ID:

  Numeric, unique identifier for each river (1 to 20).

- O2:

  Numeric, dissolved oxygen content (mg/L), higher values are better
  (positive indicator).

- PH:

  Numeric, pH value, values closer to 7 are optimal (centered
  indicator).

- germ:

  Numeric, total bacteria count, lower values are better (negative
  indicator).

- nutrient:

  Numeric, plant nutrient content (mg/L), optimal range is 10-20
  (interval indicator).

## Source

Simulated data for water quality evaluation, created for demonstration
purposes.

## Details

Water Quality Dataset

## Examples

``` r
# Load the dataset
data(water_quality)

# Preview the data
head(water_quality)
#> # A tibble: 6 × 5
#>      ID    O2    PH  germ nutrient
#>   <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1     1  4.69  6.59    51    11.9 
#> 2     2  2.03  7.86    19     6.46
#> 3     3  9.11  6.31    46     8.91
#> 4     4  8.61  7.05    46    26.4 
#> 5     5  7.13  6.5     50    23.6 
#> 6     6  2.39  6.77    38    24.6 
```
