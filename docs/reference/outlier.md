# Outlier Detection and Handling

Flags outliers in the selected numeric columns with a two-sided bound
rule and, optionally, acts on them. Detection methods: Tukey's IQR fence
(`iqr`), z-score (`zscore`, `k` standard deviations around the mean),
and the Hampel filter (`mad`, `k` scaled-MAD distances around the
median). By default the data is left untouched (`to = "detect"`): the
detection detail table is returned. With `to = "na"` detected values
become `NA` (ready for
[`impute()`](https://zhjx19.github.io/mathmodels/reference/impute.md));
with `to = "cap"` they are capped at the fence values.

## Usage

``` r
outlier(
  data,
  .cols,
  .by = NULL,
  method = c("iqr", "zscore", "mad"),
  k = NULL,
  to = c("detect", "na", "cap")
)
```

## Arguments

- data:

  A data frame.

- .cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; bounds are computed within each slice.

- method:

  `"iqr"`, `"zscore"`, or `"mad"`.

- k:

  Fence multiplier. `NULL` (default) picks the conventional value per
  method: 1.5 for `iqr`, 3 for `zscore`, 3 for `mad`.

- to:

  `"detect"` (default; data untouched, detection table returned), `"na"`
  (detected values set to `NA`), or `"cap"` (capped at the bounds).

## Value

`to = "detect"`: a tibble with one row per detected outlier (`variable`,
`.by` identifiers, `row`, `value`, `bound_low`, `bound_high`), with
attribute `summary` (per-column counts). `to = "na"` / `"cap"`: the data
frame with attribute `outliers` (the same detail table) and `summary`.

## Details

Constant columns (zero IQR / sd / MAD) cannot be scored and are skipped
with a warning.

## Examples

``` r
d = data.frame(g = rep(c("a", "b"), each = 20),
               x = c(rnorm(20), rnorm(19, 10), 100))
outlier(d, .cols = x, .by = g)
#> # A tibble: 1 × 6
#>   g     variable   row value bound_low bound_high
#> * <chr> <chr>    <int> <dbl>     <dbl>      <dbl>
#> 1 b     x           40   100      7.88       13.0
outlier(d, .cols = x, .by = g, to = "na")
#>    g           x
#> 1  a  0.45018710
#> 2  a -0.01855983
#> 3  a -0.31806837
#> 4  a -0.92936215
#> 5  a -1.48746031
#> 6  a -1.07519230
#> 7  a  1.00002880
#> 8  a -0.62126669
#> 9  a -1.38442685
#> 10 a  1.86929062
#> 11 a  0.42510038
#> 12 a -0.23864710
#> 13 a  1.05848305
#> 14 a  0.88642265
#> 15 a -0.61924305
#> 16 a  2.20610246
#> 17 a -0.25502703
#> 18 a -1.42449465
#> 19 a -0.14439960
#> 20 a  0.20753834
#> 21 b 12.30797840
#> 22 b 10.10580237
#> 23 b 10.45699881
#> 24 b  9.92284706
#> 25 b  9.66599916
#> 26 b  9.96527397
#> 27 b 10.78763961
#> 28 b 12.07524501
#> 29 b 11.02739244
#> 30 b 11.20790840
#> 31 b  8.76867658
#> 32 b 10.98389557
#> 33 b 10.21992480
#> 34 b  8.53274997
#> 35 b 10.52102274
#> 36 b  9.84124540
#> 37 b 11.46458731
#> 38 b  9.23391800
#> 39 b  9.56978825
#> 40 b          NA
outlier(d, .cols = x, .by = g, to = "cap")
#>    g           x
#> 1  a  0.45018710
#> 2  a -0.01855983
#> 3  a -0.31806837
#> 4  a -0.92936215
#> 5  a -1.48746031
#> 6  a -1.07519230
#> 7  a  1.00002880
#> 8  a -0.62126669
#> 9  a -1.38442685
#> 10 a  1.86929062
#> 11 a  0.42510038
#> 12 a -0.23864710
#> 13 a  1.05848305
#> 14 a  0.88642265
#> 15 a -0.61924305
#> 16 a  2.20610246
#> 17 a -0.25502703
#> 18 a -1.42449465
#> 19 a -0.14439960
#> 20 a  0.20753834
#> 21 b 12.30797840
#> 22 b 10.10580237
#> 23 b 10.45699881
#> 24 b  9.92284706
#> 25 b  9.66599916
#> 26 b  9.96527397
#> 27 b 10.78763961
#> 28 b 12.07524501
#> 29 b 11.02739244
#> 30 b 11.20790840
#> 31 b  8.76867658
#> 32 b 10.98389557
#> 33 b 10.21992480
#> 34 b  8.53274997
#> 35 b 10.52102274
#> 36 b  9.84124540
#> 37 b 11.46458731
#> 38 b  9.23391800
#> 39 b  9.56978825
#> 40 b 12.98515282
```
