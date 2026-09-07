# Hierarchical Clustering

Standardizes selected numeric columns (optional), builds a distance
matrix, runs hierarchical clustering
([`stats::hclust`](https://rdrr.io/r/stats/hclust.html)), and cuts the
tree into `k` clusters.

## Usage

``` r
mv_hclust(
  data,
  cols,
  k = 3,
  dist_method = "euclidean",
  method = "ward.D2",
  scale = TRUE
)
```

## Arguments

- data:

  A data frame.

- cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns to cluster on.

- k:

  Number of clusters to extract by cutting the tree.

- dist_method:

  Distance method passed to
  [`stats::dist`](https://rdrr.io/r/stats/dist.html), e.g. `"euclidean"`
  (default), `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"`, or
  `"minkowski"`.

- method:

  Agglomeration method passed to
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html), e.g.
  `"ward.D2"` (default), `"ward.D"`, `"single"`, `"complete"`,
  `"average"`, `"mcquitty"`, `"median"`, or `"centroid"`.

- scale:

  Logical. If `TRUE` (default), columns are standardized first.

## Value

An object of class `mv_hclust`: a list with components `cluster` (tibble
with `.row` and `cluster`), `heights` (tibble with `step`, `join_1`,
`join_2`, `height`; `join_1`/`join_2` use `hclust` merge encoding —
negative values index observation leaves, positive values index earlier
merged clusters), `hclust` (the raw `hclust` object, useful for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)), `data_mat`
(the standardized analysis matrix), and `meta`.

## Examples

``` r
r = mv_hclust(mtcars, cyl:carb, k = 3)
table(r$cluster$cluster)
#> 
#>  1  2  3 
#>  5 15 12 
plot_mv_dendrogram(r)
```
