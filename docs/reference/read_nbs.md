# Read and Combine National Bureau of Statistics XLS Files

This function reads multiple XLS files downloaded from the National
Bureau of Statistics of China website (<http://www.stats.gov.cn>)
without requiring renaming. It extracts the variable name from the
"indicator: XXX" text in cell A1 of each file, skips the first 3 rows of
headers, reads the first 31 rows of data, and combines the data into a
single data frame by stacking vertically based on the indicator names.
It also simplifies region names by removing suffixes such as "city",
"province", "autonomous region", "clan", or "Uyghur".

## Usage

``` r
read_nbs(paths)
```

## Arguments

- paths:

  A character vector containing the file paths to the XLS files.

## Value

A tibble containing the combined data with an additional column
`indicator` representing the original indicator names and a `region`
column with simplified region names.

## Examples

``` r
if (FALSE) { # \dontrun{
paths = c("file1.xls", "file2.xls")
data = read_nbs(paths)
} # }
```
