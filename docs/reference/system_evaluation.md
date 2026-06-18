# System Evaluation Functions for Coupling and Obstacle Analysis

These functions provide tools for system-level evaluation in
multi-indicator systems:

- `coupling_degree()`: Computes coupling degree, coordination index, and
  coupling coordination degree for subsystems.

- `obstacle_degree()`: Computes obstacle degrees for secondary
  indicators to identify key constraints in the system, enabling batch
  processing with tidyverse for grouping and summarization.

## Usage

``` r
coupling_degree(data, w = NULL, id_cols = NULL, type = "standard")

obstacle_degree(data, w = NULL, id_cols = NULL, scaled = FALSE)
```

## Arguments

- data:

  A data frame with normalized scores (usually in `[0, 1]`) as columns.

- w:

  Optional vector of weights for indicators or subsystems; defaults to
  equal weights if NULL.

- id_cols:

  Optional character vector of column names in `data` to preserve as
  identifiers (not used in calculations).

- type:

  Either "standard" for the standard coupling formula (results
  concentrated near 1) or "adjusted" for the revised formula (results
  more uniformly distributed in `[0, 1]`), as proposed by Wang Shujia,
  Kong Wei, et al. in "Misconceptions and Corrections of Domestic
  Coupling Coordination Degree Models, Journal of Natural Resources,
  2021, 36(3): 793–810 (In Chinese)"

- scaled:

  Logical. Whether to perform row normalization on obstacle degrees
  (default FALSE).

## Value

A tibble depending on the function:

- coupling_degree:

  A tibble with columns:

  - `ID`: Identifier columns specified by `id_cols` (if provided).

  - `coupling`: Coupling Degree (range 0-1).

  - `coord`: Coordination Index (range 0-1).

  - `coupling_coord`: Coupling Coordination Degree (range 0-1).

- obstacle_degree:

  A tibble with:

  - `ID`: Identifier columns specified by `id_cols` (if provided).

  - Columns for secondary indicator obstacle degrees
    (`O_{ij} = (1 - X_{ij}) * w_{ij}`).

  Suitable for grouping and summarizing (e.g., with tidyverse) to
  compute primary indicator obstacle degrees (\\ U_i \\).

## Examples

``` r
# Sample normalized subsystem scores
df = data.frame(
  ID = LETTERS[1:6],
  s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
  s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573),
  s3 = c(0.2655, 0.3721, 0.5729, 0.9082, 0.2017, 0.8984)
)
# Coupling Degree Analysis
coupling_degree(df, id_cols = "ID")        # Equal weights
#> # A tibble: 6 × 4
#>   ID    coupling coord coupling_coord
#>   <chr>    <dbl> <dbl>          <dbl>
#> 1 A        0.571 0.185          0.325
#> 2 B        0.867 0.411          0.597
#> 3 C        0.674 0.394          0.515
#> 4 D        0.967 0.687          0.815
#> 5 E        0.418 0.370          0.393
#> 6 F        0.971 0.814          0.889
coupling_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID",
                type = "adjusted")         # "adjusted" coupling degree
#> # A tibble: 6 × 4
#>   ID    coupling coord coupling_coord
#>   <chr>    <dbl> <dbl>          <dbl>
#> 1 A        0.447 0.168          0.274
#> 2 B        0.501 0.388          0.440
#> 3 C        0.455 0.409          0.432
#> 4 D        0.669 0.686          0.678
#> 5 E        0.175 0.336          0.242
#> 6 F        0.715 0.832          0.771
# Obstacle Degree Analysis
obstacle_degree(df, id_cols = "ID")        # Equal weights
#> # A tibble: 6 × 4
#>   ID         s1     s2     s3
#>   <chr>   <dbl>  <dbl>  <dbl>
#> 1 A     0.328   0.243  0.245 
#> 2 B     0.274   0.106  0.209 
#> 3 C     0.150   0.314  0.142 
#> 4 D     0.109   0.173  0.0306
#> 5 E     0.326   0.0370 0.266 
#> 6 F     0.00417 0.148  0.0339
obstacle_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID")
#> # A tibble: 6 × 4
#>   ID         s1     s2     s3
#>   <chr>   <dbl>  <dbl>  <dbl>
#> 1 A     0.394   0.218  0.220 
#> 2 B     0.329   0.0953 0.188 
#> 3 C     0.180   0.282  0.128 
#> 4 D     0.131   0.156  0.0275
#> 5 E     0.392   0.0333 0.239 
#> 6 F     0.00500 0.133  0.0305
```
