# Clustering Validity Metrics

Computes validity metrics for a partition of the selected numeric
columns: per-observation silhouette widths plus the average silhouette,
the Davies-Bouldin index, the Calinski-Harabasz index, and the Dunn
index. All metrics are computed in base R from the (optionally
standardized) analysis matrix and its distances.

## Usage

``` r
mv_cluster_metrics(
  data,
  cols,
  cluster,
  dist_method = "euclidean",
  scale = TRUE
)
```

## Arguments

- data:

  A data frame.

- cols:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric columns used for distances.

- cluster:

  A vector of cluster labels with one value per row of `data` (factor,
  character, or numeric).

- dist_method:

  Distance method passed to
  [`stats::dist`](https://rdrr.io/r/stats/dist.html).

- scale:

  Logical. If `TRUE` (default), columns are standardized first.

## Value

An object of class `mv_cluster_metrics`: a list with components
`silhouette` (tibble with `.row`, `cluster`, `neighbor`, `sil_width`),
`summary` (long tibble with `metric` and `value`: `avg_silhouette`,
`davies_bouldin`, `calinski_harabasz`, `dunn`), `data_mat` (the
standardized analysis matrix), and `meta`.

## Examples

``` r
set.seed(1)
km = mv_kmeans(iris, Sepal.Length:Petal.Width, centers = 3)
m = mv_cluster_metrics(iris, Sepal.Length:Petal.Width, km$cluster$cluster)
m$summary
#> # A tibble: 4 × 2
#>   metric               value
#>   <chr>                <dbl>
#> 1 avg_silhouette      0.460 
#> 2 davies_bouldin      0.834 
#> 3 calinski_harabasz 242.    
#> 4 dunn                0.0265
plot_mv_silhouette(m)
```
