# Regional Economics Functions

A collection of functions for calculating regional economics indices,
including Location Quotient (LQ), Herfindahl-Hirschman Index (HHI), and
Ellison-Glaeser Index (EG). These functions are designed for analyzing
regional and industrial data to assess spatial concentration and market
structure.

## Usage

``` r
LQ(data, region, total, .cols, .by = NULL)

HHI(x, scaled = FALSE)

EG(data, region, industry, y)
```

## Arguments

- data:

  A data frame containing the necessary data for calculations. For `LQ`,
  the first row is assumed to be the national total, with the first two
  columns specifying the region and total, optionally including a
  grouping column. For `EG`, it contains region, industry, and indicator
  columns (e.g., employment or output).

- region:

  The name of the column in `data` specifying the region (used in `LQ`
  and `EG`).

- total:

  The name of the column in `data` specifying the total (used in `LQ`).

- .cols:

  The columns in `data` for which to calculate the location quotient
  (used in `LQ`).

- .by:

  Optional grouping column for `LQ`, defaults to `NULL` (no grouping).

- x:

  A numeric vector for calculating the HHI (used in `HHI`).

- scaled:

  Logical; if `TRUE`, the HHI is scaled to account for the number of
  firms. Defaults to `FALSE`.

- industry:

  The name of the column in `data` specifying the industry (used in
  `EG`).

- y:

  The name of the column in `data` specifying the indicator (e.g.,
  employment or output, used in `EG`).

## Value

- LQ:

  A data frame with the region column and calculated location quotients
  for the specified columns, optionally grouped by `.by`.

- HHI:

  A numeric value representing the HHI, either scaled or unscaled.

- EG:

  A tibble with two columns: the industry name and the corresponding EG
  index value.

## Details

- LQ:

  Calculates the Location Quotient for multiple columns with optional
  grouping. The LQ measures the relative concentration of an industry in
  a region compared to a national benchmark. The function assumes the
  first row of `data` contains national totals, with the first two
  columns specifying the region and total, and the remaining columns
  used for LQ calculation.

- HHI:

  Calculates the Herfindahl-Hirschman Index, a measure of market
  concentration based on the squared sum of market shares. If
  `scaled = TRUE`, the HHI is normalized to account for the number of
  firms.

- EG:

  Calculates the Ellison-Glaeser Index, which measures the geographic
  concentration of an industry while controlling for firm size
  distribution and random distribution effects. It requires data on
  regions, industries, and an indicator (e.g., employment or output).

## Examples

``` r
# Example data
data = data.frame(
  region = c("National", "Region_A", "Region_B"),
  total = c(10000, 4000, 6000),
  industry1 = c(2000, 1000, 1000),
  industry2 = c(6000, 2000, 4000)
)

# Calculate Location Quotient
LQ(data, region, total, starts_with("industry"))
#>     region industry1 industry2
#> 1 Region_A 1.2500000 0.8333333
#> 2 Region_B 0.8333333 1.1111111

# Calculate HHI
x = c(50, 30, 20)
# Calculate the raw HHI
HHI(x)
#> [1] 0.38
# Calculate the standard (scaled) HHI
HHI(x, scaled = TRUE)
#> [1] 0.07

# Example data for EG
eg_data = data.frame(
 region = c("R1", "R1", "R1", "R1", "R2", "R2", "R3", "R3", "R1", "R2", "R2", "R3"),
 industry = c("A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B"),
 employment = c(250, 200, 150, 100, 20, 15, 10, 5, 50, 200, 150, 50)
)
EG(eg_data, region, industry, y = employment)
#> # A tibble: 2 × 2
#>   industry    EG
#>   <chr>    <dbl>
#> 1 A        0.131
#> 2 B        0.918
```
