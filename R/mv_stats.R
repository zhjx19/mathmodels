# =============================================================================
# Multivariate statistics toolkit for mathmodels
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Validate that the input is a data frame.
.mv_frame = function(data) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  invisible(TRUE)
}

# Build the numeric analysis matrix from an already tidy-selected data frame,
# checking column types and dropping NA rows with bookkeeping. Inputs without
# meaningful row names (e.g. tibbles) get automatic numeric row labels so that
# downstream outputs (dendrograms, distance-based plots) stay labelled.
.mv_matrix = function(sub, arg = "cols", min_cols = 1L) {
  nm = names(sub)
  if (length(nm) < min_cols) {
    stop(arg, " must select at least ", min_cols, " column(s).", call. = FALSE)
  }
  bad = nm[!vapply(sub, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(
      "Selected columns must be numeric. Non-numeric: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  keep = stats::complete.cases(sub)
  # Inputs without meaningful row names (tibbles, auto-named) get numeric
  # labels; convert before subsetting because tibble `[` drops row names.
  rn = rownames(sub)
  auto_named = is.null(rn) || all(rn == "") ||
    identical(attr(sub, "rownames"), paste0("auto:", seq_len(nrow(sub)))) ||
    identical(as.character(rn), as.character(seq_len(nrow(sub))))
  mat = as.matrix(sub)[keep, , drop = FALSE]
  if (auto_named) {
    rownames(mat) = as.character(which(keep))
  }
  list(
    mat = mat,
    names = nm,
    row = which(keep),
    n_used = sum(keep),
    n_dropped = sum(!keep)
  )
}

.mv_check_constant = function(mat) {
  s = apply(mat, 2, stats::sd)
  bad = names(s)[s == 0 | is.na(s)]
  if (length(bad) > 0) {
    stop(
      "Cannot standardize constant column(s): ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.mv_scale = function(mat) {
  scale(mat, center = TRUE, scale = TRUE)
}

.mv_dist_mat = function(mat, method = "euclidean") {
  as.matrix(stats::dist(mat, method = method))
}

# Rousseeuw silhouette widths from a full distance matrix and integer codes.
.mv_silhouette = function(d, codes) {
  n = nrow(d)
  uni = sort(unique(codes))
  s = numeric(n)
  nb = character(n)
  for (i in seq_len(n)) {
    ci = codes[i]
    own = codes == ci
    own[i] = FALSE
    if (!any(own)) {
      s[i] = 0
      b = Inf
      best = ci
    } else {
      a = mean(d[i, own])
      b = Inf
      best = ci
      for (cj in uni[uni != ci]) {
        m = mean(d[i, codes == cj])
        if (m < b) {
          b = m
          best = cj
        }
      }
      denom = max(a, b)
      s[i] = if (denom > 0) (b - a) / denom else 0
    }
    nb[i] = as.character(best)
  }
  list(width = s, neighbor = nb, avg = mean(s))
}

# Silhouette as a tidy tibble from an analysis matrix and label vector.
# Neighbor codes are mapped back to the original label space.
.mv_sil_df = function(mat, labels, dist_method = "euclidean") {
  f = factor(labels)
  codes = as.integer(f)
  d = .mv_dist_mat(mat, dist_method)
  sil = .mv_silhouette(d, codes)
  tibble::tibble(
    .row = seq_along(codes),
    cluster = f,
    neighbor = levels(f)[as.integer(sil$neighbor)],
    sil_width = sil$width
  )
}

# Within-cluster sum of squares from an analysis matrix and integer codes.
.mv_wss = function(mat, codes) {
  w = 0
  for (g in unique(codes)) {
    sub = mat[codes == g, , drop = FALSE]
    if (nrow(sub) > 0) {
      w = w + sum(scale(sub, center = TRUE, scale = FALSE)^2)
    }
  }
  w
}

# Davies-Bouldin index (euclidean centroids in the analysis space); cluster
# scatter is the mean euclidean distance of members to the cluster centroid,
# following Davies & Bouldin (1979).
.mv_db = function(mat, codes) {
  uni = sort(unique(codes))
  k = length(uni)
  p = ncol(mat)
  cen = t(vapply(uni, function(g) {
    colMeans(mat[codes == g, , drop = FALSE])
  }, numeric(p)))
  scatter = vapply(seq_len(k), function(i) {
    xg = mat[codes == uni[i], , drop = FALSE]
    mean(sqrt(rowSums((xg - matrix(cen[i, ], nrow(xg), ncol(xg), byrow = TRUE))^2)))
  }, numeric(1))
  r = matrix(0, k, k)
  for (i in seq_len(k - 1)) {
    for (j in seq.int(i + 1, k)) {
      dd = sqrt(sum((cen[i, ] - cen[j, ])^2))
      r[i, j] = r[j, i] = (scatter[i] + scatter[j]) / dd
    }
  }
  mean(vapply(seq_len(k), function(i) max(r[i, -i]), numeric(1)))
}

# Calinski-Harabasz index (euclidean in the analysis space).
.mv_ch = function(mat, codes) {
  n = nrow(mat)
  uni = sort(unique(codes))
  k = length(uni)
  p = ncol(mat)
  grand = colMeans(mat)
  cen = t(vapply(uni, function(g) {
    colMeans(mat[codes == g, , drop = FALSE])
  }, numeric(p)))
  w = .mv_wss(mat, codes)
  b = sum(vapply(seq_len(k), function(i) {
    sum(codes == uni[i]) * sum((cen[i, ] - grand)^2)
  }, numeric(1)))
  (b / (k - 1)) / (w / (n - k))
}

# Dunn index: min inter-cluster distance / max intra-cluster distance.
.mv_dunn = function(d, codes) {
  uni = sort(unique(codes))
  k = length(uni)
  inter = Inf
  for (i in seq_len(k - 1)) {
    for (j in seq.int(i + 1, k)) {
      sub = d[codes == uni[i], codes == uni[j], drop = FALSE]
      inter = min(inter, min(sub))
    }
  }
  intra = 0
  for (g in uni) {
    sub = d[codes == g, codes == g, drop = FALSE]
    if (nrow(sub) > 1) intra = max(intra, max(sub))
  }
  if (intra == 0) Inf else inter / intra
}

# Segments and labels for drawing a dendrogram with ggplot2.
.mv_dendro_segments = function(hc) {
  n = length(hc$order)
  merge = hc$merge
  h = hc$height
  xleaf = integer(n)
  xleaf[hc$order] = seq_len(n)
  xpos = numeric(n - 1)
  getx = function(v) if (v < 0) xleaf[-v] else xpos[v]
  gety = function(v) if (v < 0) 0 else h[v]
  segs = vector("list", 3 * (n - 1))
  m = 1
  for (i in seq_len(n - 1)) {
    a = merge[i, 1]
    b = merge[i, 2]
    xa = getx(a)
    xb = getx(b)
    xpos[i] = (xa + xb) / 2
    segs[[m]] = data.frame(x = xa, y = gety(a), xend = xa, yend = h[i]); m = m + 1
    segs[[m]] = data.frame(x = xb, y = gety(b), xend = xb, yend = h[i]); m = m + 1
    segs[[m]] = data.frame(x = xa, y = h[i], xend = xb, yend = h[i]); m = m + 1
  }
  segments = do.call(rbind, segs)
  labels = tibble::tibble(x = seq_len(n), label = hc$labels[hc$order])
  list(segments = segments, labels = labels)
}

# -----------------------------------------------------------------------------
# Principal component analysis
# -----------------------------------------------------------------------------

#' Principal Component Analysis
#'
#' Performs a PCA on selected numeric columns and returns tidy components:
#' observation scores, variable loadings, and the variance decomposition.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to analyse.
#' @param scale Logical. If `TRUE` (default), columns are standardized to unit
#'   variance before decomposition; recommended when indicators use mixed units.
#' @param center Logical. If `TRUE` (default), columns are centered.
#'
#' @return An object of class `mv_pca`: a list with components
#'   `scores` (tibble with `.row` and `PC1`...), `loadings` (tibble with
#'   `variable` and `PC1`...), `variance` (tibble with `PC`, `sdev`,
#'   `variance`, `prop_var`, `cum_var`), and `meta`.
#'
#' @examples
#' r = mv_pca(mtcars, cyl:carb)
#' r$variance
#' plot_mv_scree(r)
#' @export
mv_pca = function(data, cols, scale = TRUE, center = TRUE) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  mat = prep$mat
  if (!is.logical(scale) || length(scale) != 1 || !is.logical(center) ||
      length(center) != 1) {
    stop("scale and center must be single logical values.", call. = FALSE)
  }
  if (scale) .mv_check_constant(mat)
  fit = stats::prcomp(mat, center = center, scale. = scale, retx = TRUE)
  k = length(fit$sdev)
  scores = tibble::tibble(.row = prep$row, tibble::as_tibble(fit$x))
  loadings = tibble::tibble(variable = prep$names, tibble::as_tibble(fit$rotation))
  variance = tibble::tibble(
    PC = paste0("PC", seq_len(k)),
    sdev = fit$sdev,
    variance = fit$sdev^2,
    prop_var = fit$sdev^2 / sum(fit$sdev^2),
    cum_var = cumsum(fit$sdev^2 / sum(fit$sdev^2))
  )
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    n_pc = k,
    center = center,
    scale = scale,
    variables = list(prep$names)
  )
  structure(
    list(scores = scores, loadings = loadings, variance = variance, meta = meta),
    class = "mv_pca"
  )
}

#' @export
print.mv_pca = function(x, ...) {
  cat("Principal component analysis (mv_pca)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA)\n\n",
    x$meta$n_used, x$meta$n_dropped
  ))
  print(x$variance)
  cat("\nLoadings:\n")
  print(x$loadings)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Exploratory factor analysis
# -----------------------------------------------------------------------------

#' Exploratory Factor Analysis
#'
#' Performs an exploratory factor analysis on selected numeric columns via
#' maximum likelihood (`stats::factanal`) on the standardized data, and
#' returns tidy loadings, uniquenesses, variance shares, and (optionally)
#' factor scores.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to analyse.
#' @param factors Number of factors to extract.
#' @param rotation Character. Rotation passed to `factanal`: one of
#'   `"varimax"` (default), `"promax"`, `"none"`, or the name of a rotation
#'   function.
#' @param scores Either `"regression"` (default) to compute factor scores, or
#'   `"none"`.
#'
#' @return An object of class `mv_efa`: a list with components `loadings`
#'   (tibble with `variable`, `Factor1`..., and `uniqueness`), `uniquenesses`,
#'   `variance` (tibble with `factor`, `ss_loadings`, `prop_var`, `cum_var`),
#'   `scores` (tibble with `.row` and `Factor1`..., or `NULL`), and `meta`.
#'
#' @examples
#' r = mv_efa(mtcars, cyl:carb, factors = 2)
#' r$loadings
#' @export
mv_efa = function(data, cols, factors, rotation = "varimax",
                  scores = c("regression", "none")) {
  scores = match.arg(scores)
  if (!is.character(rotation) || length(rotation) != 1) {
    stop("rotation must be a single character value.", call. = FALSE)
  }
  known_rot = c("varimax", "promax", "none")
  if (!(rotation %in% known_rot) &&
      !exists(rotation, mode = "function", envir = .GlobalEnv)) {
    stop(
      "rotation must be one of 'varimax', 'promax', 'none', ",
      "or the name of an available rotation function.",
      call. = FALSE
    )
  }
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  p = length(prep$names)
  if (!is.numeric(factors) || length(factors) != 1 || factors < 1 ||
      factors != floor(factors)) {
    stop("factors must be a positive integer.", call. = FALSE)
  }
  if (factors >= p) {
    stop("factors must be smaller than the number of selected variables.",
         call. = FALSE)
  }
  if (prep$n_used <= p) {
    stop("Need more complete observations than selected variables.",
         call. = FALSE)
  }
  .mv_check_constant(prep$mat)
  mat = .mv_scale(prep$mat)
  fit = stats::factanal(
    mat,
    factors = factors,
    rotation = rotation,
    scores = scores
  )
  ld = unclass(fit$loadings)
  f_names = colnames(ld)
  loadings = tibble::tibble(
    variable = rownames(ld),
    tibble::as_tibble(ld),
    uniqueness = as.numeric(fit$uniquenesses)
  )
  uniquenesses = tibble::tibble(
    variable = rownames(ld),
    uniqueness = as.numeric(fit$uniquenesses)
  )
  ss = colSums(ld^2)
  variance = tibble::tibble(
    factor = f_names,
    ss_loadings = ss,
    prop_var = ss / p,
    cum_var = cumsum(ss / p)
  )
  f_scores = if (!is.null(fit$scores)) {
    tibble::tibble(.row = prep$row, tibble::as_tibble(fit$scores))
  } else {
    NULL
  }
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    factors = factors,
    rotation = rotation,
    variables = list(prep$names)
  )
  structure(
    list(
      loadings = loadings,
      uniquenesses = uniquenesses,
      variance = variance,
      scores = f_scores,
      meta = meta
    ),
    class = "mv_efa"
  )
}

#' @export
print.mv_efa = function(x, ...) {
  cat("Exploratory factor analysis (mv_efa)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); factors: %d; rotation: %s\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$factors, x$meta$rotation
  ))
  print(x$variance)
  cat("\nLoadings:\n")
  print(x$loadings)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Hierarchical clustering
# -----------------------------------------------------------------------------

#' Hierarchical Clustering
#'
#' Standardizes selected numeric columns (optional), builds a distance matrix,
#' runs hierarchical clustering (`stats::hclust`), and cuts the tree into `k`
#' clusters.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to cluster on.
#' @param k Number of clusters to extract by cutting the tree.
#' @param dist_method Distance method passed to `stats::dist`, e.g.
#'   `"euclidean"` (default), `"maximum"`, `"manhattan"`, `"canberra"`,
#'   `"binary"`, or `"minkowski"`.
#' @param method Agglomeration method passed to `stats::hclust`, e.g.
#'   `"ward.D2"` (default), `"ward.D"`, `"single"`, `"complete"`, `"average"`,
#'   `"mcquitty"`, `"median"`, or `"centroid"`.
#' @param scale Logical. If `TRUE` (default), columns are standardized first.
#'
#' @return An object of class `mv_hclust`: a list with components `cluster`
#'   (tibble with `.row` and `cluster`), `heights` (tibble with `step`,
#'   `join_1`, `join_2`, `height`; `join_1`/`join_2` use `hclust` merge
#'   encoding — negative values index observation leaves, positive values
#'   index earlier merged clusters), `hclust` (the raw `hclust` object,
#'   useful for `plot()`), `data_mat` (the standardized analysis matrix), and
#'   `meta`.
#'
#' @examples
#' r = mv_hclust(mtcars, cyl:carb, k = 3)
#' table(r$cluster$cluster)
#' plot_mv_dendrogram(r)
#' @export
mv_hclust = function(data, cols, k = 3, dist_method = "euclidean",
                     method = "ward.D2", scale = TRUE) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  dist_method = match.arg(dist_method, c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  ))
  method = match.arg(method, c(
    "ward.D", "ward.D2", "single", "complete", "average",
    "mcquitty", "median", "centroid"
  ))
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  if (k >= prep$n_used) {
    stop("k must be smaller than the number of complete observations.",
         call. = FALSE)
  }
  if (scale) .mv_check_constant(prep$mat)
  mat = if (scale) .mv_scale(prep$mat) else prep$mat
  d = stats::dist(mat, method = dist_method)
  hc = stats::hclust(d, method = method)
  cl = stats::cutree(hc, k = k)
  n_step = nrow(hc$merge)
  heights = tibble::tibble(
    step = seq_len(n_step),
    join_1 = hc$merge[, 1],
    join_2 = hc$merge[, 2],
    height = hc$height
  )
  cluster = tibble::tibble(.row = prep$row, cluster = factor(cl))
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    k = k,
    dist_method = dist_method,
    method = method,
    scale = scale
  )
  structure(
    list(
      cluster = cluster,
      heights = heights,
      hclust = hc,
      data_mat = mat,
      meta = meta
    ),
    class = "mv_hclust"
  )
}

#' @export
print.mv_hclust = function(x, ...) {
  cat("Hierarchical clustering (mv_hclust)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); k = %d; method = %s; distance = %s\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$k, x$meta$method,
    x$meta$dist_method
  ))
  print(table(x$cluster$cluster))
  invisible(x)
}

# -----------------------------------------------------------------------------
# K-means clustering
# -----------------------------------------------------------------------------

#' K-means Clustering
#'
#' Standardizes selected numeric columns (optional) and runs `stats::kmeans`,
#' returning tidy cluster assignments, centers, and within-cluster sums of
#' squares.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to cluster on.
#' @param centers Either the number of clusters (`k`, a single integer) or a
#'   numeric matrix of initial cluster centers with the same number of
#'   columns as `cols`, passed to `stats::kmeans`. When a named matrix of
#'   centers is supplied, its row names are kept as cluster labels.
#' @param nstart Number of random starts passed to `stats::kmeans`; ignored
#'   when `centers` is a matrix.
#' @param scale Logical. If `TRUE` (default), columns are standardized first.
#' @param iter.max Maximum iterations passed to `stats::kmeans`.
#'
#' @return An object of class `mv_kmeans`: a list with components `cluster`
#'   (tibble with `.row` and `cluster`), `centers` (tibble with `cluster` and
#'   one column per variable), `within_ss` (tibble with `cluster`, `size`,
#'   `withinss`, `share`), `data_mat` (the standardized analysis matrix), and
#'   `meta`.
#'
#' @examples
#' r = mv_kmeans(iris, Sepal.Length:Petal.Width, centers = 3)
#' r$within_ss
#' plot_mv_silhouette(r)
#' @export
mv_kmeans = function(data, cols, centers, nstart = 25, scale = TRUE,
                     iter.max = 100) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  if (!is.numeric(nstart) || length(nstart) != 1 || nstart < 1) {
    stop("nstart must be a positive integer.", call. = FALSE)
  }
  center_names = NULL
  if (is.matrix(centers)) {
    if (!is.numeric(centers)) {
      stop("centers must be a number of clusters or a numeric matrix.",
           call. = FALSE)
    }
    if (ncol(centers) != ncol(prep$mat)) {
      stop(
        "centers must have one column per selected variable (",
        ncol(prep$mat), ").",
        call. = FALSE
      )
    }
    if (is.null(rownames(centers)) || all(rownames(centers) == "")) {
      rownames(centers) = paste0("center_", seq_len(nrow(centers)))
    }
    center_names = rownames(centers)
  } else {
    if (!is.numeric(centers) || length(centers) != 1 || centers < 1 ||
        centers != floor(centers)) {
      stop(
        "centers must be a positive integer (number of clusters) or a matrix.",
        call. = FALSE
      )
    }
  }
  if (scale) .mv_check_constant(prep$mat)
  mat = if (scale) .mv_scale(prep$mat) else prep$mat
  # With a matrix of centers, kmeans treats nstart > 1 as multiple sets of
  # initial centers; force nstart = 1 so the user's centers are honoured.
  fit = stats::kmeans(
    mat,
    centers = centers,
    nstart = if (is.matrix(centers)) 1L else nstart,
    iter.max = iter.max
  )
  k = nrow(fit$centers)
  if (is.null(center_names)) {
    center_names = as.character(seq_len(k))
  } else {
    center_names = center_names[seq_len(k)]
  }
  center_f = factor(center_names, levels = center_names)
  cluster_levels = center_f[fit$cluster]
  cluster = tibble::tibble(.row = prep$row, cluster = cluster_levels)
  centers_tb = tibble::tibble(
    cluster = factor(center_names, levels = center_names),
    tibble::as_tibble(fit$centers)
  )
  within_ss = tibble::tibble(
    cluster = factor(center_names, levels = center_names),
    size = as.integer(fit$size),
    withinss = fit$withinss,
    share = fit$withinss / fit$tot.withinss
  )
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    k = k,
    totss = fit$totss,
    tot_withinss = fit$tot.withinss,
    betweenss = fit$betweenss,
    iter = fit$iter,
    scale = scale
  )
  structure(
    list(
      cluster = cluster,
      centers = centers_tb,
      within_ss = within_ss,
      data_mat = mat,
      meta = meta
    ),
    class = "mv_kmeans"
  )
}

#' @export
print.mv_kmeans = function(x, ...) {
  cat("K-means clustering (mv_kmeans)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); k = %d\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$k
  ))
  print(table(x$cluster$cluster))
  cat(sprintf(
    "\nTotal within SS: %.3f; between SS: %.3f; total SS: %.3f\n",
    x$meta$tot_withinss, x$meta$betweenss, x$meta$totss
  ))
  invisible(x)
}

# -----------------------------------------------------------------------------
# Clustering validity metrics
# -----------------------------------------------------------------------------

#' Clustering Validity Metrics
#'
#' Computes validity metrics for a partition of the selected numeric columns:
#' per-observation silhouette widths plus the average silhouette, the
#' Davies-Bouldin index, the Calinski-Harabasz index, and the Dunn index.
#' All metrics are computed in base R from the (optionally standardized)
#' analysis matrix and its distances.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns used for distances.
#' @param cluster A vector of cluster labels with one value per row of `data`
#'   (factor, character, or numeric).
#' @param dist_method Distance method passed to `stats::dist`.
#' @param scale Logical. If `TRUE` (default), columns are standardized first.
#'
#' @return An object of class `mv_cluster_metrics`: a list with components
#'   `silhouette` (tibble with `.row`, `cluster`, `neighbor`, `sil_width`),
#'   `summary` (long tibble with `metric` and `value`:
#'   `avg_silhouette`, `davies_bouldin`, `calinski_harabasz`, `dunn`),
#'   `data_mat` (the standardized analysis matrix), and `meta`.
#'
#' @examples
#' set.seed(1)
#' km = mv_kmeans(iris, Sepal.Length:Petal.Width, centers = 3)
#' m = mv_cluster_metrics(iris, Sepal.Length:Petal.Width, km$cluster$cluster)
#' m$summary
#' plot_mv_silhouette(m)
#' @export
mv_cluster_metrics = function(data, cols, cluster, dist_method = "euclidean",
                              scale = TRUE) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  dist_method = match.arg(dist_method, c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  ))
  if (length(cluster) != nrow(data)) {
    stop("cluster must be a vector with one label per row of data.",
         call. = FALSE)
  }
  labels = cluster[prep$row]
  if (anyNA(labels)) {
    stop("cluster must not contain NA in the analysed rows.", call. = FALSE)
  }
  f = factor(labels)
  k = nlevels(f)
  if (k < 2) {
    stop("cluster must contain at least 2 distinct labels.", call. = FALSE)
  }
  if (k >= prep$n_used) {
    stop("cluster must contain fewer distinct labels than observations.",
         call. = FALSE)
  }
  if (scale) .mv_check_constant(prep$mat)
  mat = if (scale) .mv_scale(prep$mat) else prep$mat
  codes = as.integer(f)
  d = .mv_dist_mat(mat, dist_method)
  sil = .mv_silhouette(d, codes)
  silhouette = tibble::tibble(
    .row = prep$row,
    cluster = f,
    neighbor = levels(f)[as.integer(sil$neighbor)],
    sil_width = sil$width
  )
  summary = tibble::tibble(
    metric = c(
      "avg_silhouette", "davies_bouldin", "calinski_harabasz", "dunn"
    ),
    value = c(sil$avg, .mv_db(mat, codes), .mv_ch(mat, codes),
              .mv_dunn(d, codes))
  )
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    k = k,
    dist_method = dist_method,
    scale = scale
  )
  structure(
    list(
      silhouette = silhouette,
      summary = summary,
      data_mat = mat,
      meta = meta
    ),
    class = "mv_cluster_metrics"
  )
}

#' @export
print.mv_cluster_metrics = function(x, ...) {
  cat("Clustering validity metrics (mv_cluster_metrics)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); k = %d\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$k
  ))
  print(x$summary)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Discriminant analysis
# -----------------------------------------------------------------------------

#' Linear and Quadratic Discriminant Analysis
#'
#' Fits LDA or QDA (`MASS::lda` / `MASS::qda`) on selected numeric columns,
#' reports training predictions with a confusion matrix and accuracy, and can
#' classify new observations.
#'
#' @param data A data frame.
#' @param group <[`tidy-select`][tidyselect::language]> The single grouping
#'   column (coerced to a factor; must not contain `NA`).
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric predictor
#'   columns; must not include the grouping column (exclude it with `-group`).
#' @param method Either `"lda"` (default) or `"qda"`.
#' @param prior Optional prior probabilities for the groups, passed to
#'   `MASS::lda()` / `MASS::qda()`.
#' @param newdata Optional data frame with the same predictor columns to
#'   classify; see `new_predictions` in the return value.
#'
#' @return An object of class `mv_discrim`: a list with components
#'   `predictions` (tibble with `.row`, `actual`, `predicted`, and
#'   `post_<level>` posterior columns), `confusion` (long tibble with
#'   `actual`, `predicted`, `n`), `scores` (tibble with `.row` and
#'   `LD1`... — LDA only, `NULL` for QDA), `structure` (tibble with
#'   `variable` and `LD1`... for LDA, or per-group means for QDA),
#'   `new_predictions` (tibble with `.row` — the row number within `newdata`
#'   — plus `predicted` and posterior columns; rows of `newdata` with missing
#'   predictor values are skipped), and `meta`.
#'
#' @examples
#' r = mv_discrim(iris, Species, Sepal.Length:Petal.Width)
#' r$meta$accuracy
#' r$confusion
#' @export
mv_discrim = function(data, group, cols, method = c("lda", "qda"),
                      prior = NULL, newdata = NULL) {
  method = match.arg(method)
  .mv_frame(data)
  gsub = dplyr::select(data, {{ group }})
  if (ncol(gsub) != 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  gname = names(gsub)
  sub = dplyr::select(data, {{ cols }})
  if (gname %in% names(sub)) {
    stop(
      "The grouping column '", gname, "' must not be among the predictors; ",
      "exclude it from cols, e.g. cols = -", gname, ".",
      call. = FALSE
    )
  }
  prep = .mv_matrix(sub)
  gv = gsub[[1]]
  if (anyNA(gv)) {
    stop("group must not contain NA.", call. = FALSE)
  }
  grp = droplevels(as.factor(gv[prep$row]))
  if (nlevels(grp) < 2) {
    stop("group must have at least 2 levels in the complete rows.",
         call. = FALSE)
  }
  X = as.data.frame(prep$mat)
  fit = if (method == "lda") {
    if (is.null(prior)) {
      MASS::lda(x = X, grouping = grp)
    } else {
      MASS::lda(x = X, grouping = grp, prior = prior)
    }
  } else {
    if (is.null(prior)) {
      MASS::qda(x = X, grouping = grp)
    } else {
      MASS::qda(x = X, grouping = grp, prior = prior)
    }
  }
  ptr = stats::predict(fit, newdata = X)
  post = as.data.frame(ptr$posterior)
  names(post) = paste0("post_", names(post))
  predictions = tibble::tibble(
    .row = prep$row,
    actual = grp,
    predicted = factor(as.character(ptr$class), levels = levels(grp)),
    post
  )
  confusion = tibble::as_tibble(as.data.frame(
    table(actual = predictions$actual, predicted = predictions$predicted),
    responseName = "n"
  ))
  f_scores = NULL
  structure_tb = NULL
  if (method == "lda") {
    if (!is.null(ptr$x)) {
      f_scores = tibble::tibble(.row = prep$row, tibble::as_tibble(ptr$x))
    }
    scaling = fit$scaling
    structure_tb = tibble::tibble(
      variable = rownames(scaling),
      tibble::as_tibble(scaling)
    )
  } else {
    means = fit$means
    structure_tb = tibble::tibble(
      group = factor(rownames(means), levels = levels(grp)),
      tibble::as_tibble(means)
    )
  }
  new_pred = NULL
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) {
      stop("newdata must be a data.frame.", call. = FALSE)
    }
    missing = setdiff(prep$names, names(newdata))
    if (length(missing) > 0) {
      stop(
        "newdata is missing predictor column(s): ",
        paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
    nd = newdata[prep$names]
    bad = prep$names[!vapply(nd, is.numeric, logical(1))]
    if (length(bad) > 0) {
      stop(
        "newdata predictor columns must be numeric. Non-numeric: ",
        paste(bad, collapse = ", "), ".",
        call. = FALSE
      )
    }
    keep = stats::complete.cases(nd)
    ndm = as.matrix(nd[keep, , drop = FALSE])
    colnames(ndm) = prep$names
    pnd = stats::predict(fit, newdata = as.data.frame(ndm))
    post_nd = as.data.frame(pnd$posterior)
    names(post_nd) = paste0("post_", names(post_nd))
    new_pred = tibble::tibble(
      .row = which(keep),
      predicted = factor(as.character(pnd$class), levels = levels(grp)),
      post_nd
    )
  }
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    method = method,
    accuracy = mean(predictions$actual == predictions$predicted),
    n_predictors = length(prep$names),
    levels = list(levels(grp)),
    prior = list(fit$prior)
  )
  structure(
    list(
      predictions = predictions,
      confusion = confusion,
      scores = f_scores,
      structure = structure_tb,
      new_predictions = new_pred,
      meta = meta
    ),
    class = "mv_discrim"
  )
}

#' @export
print.mv_discrim = function(x, ...) {
  cat(sprintf("Discriminant analysis (mv_discrim, method = %s)\n", x$meta$method))
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); accuracy: %.3f\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$accuracy
  ))
  print(x$confusion)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Correspondence analysis
# -----------------------------------------------------------------------------

#' Simple Correspondence Analysis
#'
#' Builds a contingency table from two columns of a data frame and runs simple
#' correspondence analysis (`MASS::corresp`), returning row and column
#' coordinates and the inertia decomposition.
#'
#' @param data A data frame.
#' @param row <[`tidy-select`][tidyselect::language]> The single row-factor column.
#' @param col <[`tidy-select`][tidyselect::language]> The single column-factor column.
#' @param nf Number of dimensions to compute; must be smaller than both the
#'   number of rows and the number of columns of the table.
#'
#' @return An object of class `mv_corresp`: a list with components
#'   `row_coord` (tibble with `row` and `Dim1`...), `col_coord` (tibble with
#'   `col` and `Dim1`...), `inertia` (tibble with `dimension`, `inertia`,
#'   `prop`, `cum`), and `meta`.
#'
#' @examples
#' set.seed(42)
#' df = data.frame(
#'   breed = rep(c("A", "B", "C"), each = 30),
#'   grade = sample(c("low", "mid", "high"), 90, TRUE)
#' )
#' r = mv_corresp(df, breed, grade)
#' r$inertia
#' @export
mv_corresp = function(data, row, col, nf = 2) {
  .mv_frame(data)
  rsub = dplyr::select(data, {{ row }})
  csub = dplyr::select(data, {{ col }})
  if (ncol(rsub) != 1) {
    stop("row must select exactly one column.", call. = FALSE)
  }
  if (ncol(csub) != 1) {
    stop("col must select exactly one column.", call. = FALSE)
  }
  if (identical(names(rsub), names(csub))) {
    stop("row and col must be different columns.", call. = FALSE)
  }
  rv = rsub[[1]]
  cv = csub[[1]]
  keep = !is.na(rv) & !is.na(cv)
  rf = factor(rv[keep])
  cf = factor(cv[keep])
  tab = table(rf, cf)
  nr = nrow(tab)
  nc = ncol(tab)
  if (nr < 2 || nc < 2) {
    stop("Both row and col need at least 2 levels after dropping NA.",
         call. = FALSE)
  }
  if (!is.numeric(nf) || length(nf) != 1 || nf < 1 || nf != floor(nf)) {
    stop("nf must be a positive integer.", call. = FALSE)
  }
  if (nf > min(nr, nc) - 1) {
    stop("nf must be smaller than min(nrow(table), ncol(table)).",
         call. = FALSE)
  }
  fit = MASS::corresp(unclass(tab), nf = nf)
  dims = paste0("Dim", seq_len(nf))
  rc = fit$rscore
  cs = fit$cscore
  colnames(rc) = dims
  colnames(cs) = dims
  row_coord = tibble::tibble(row = rownames(rc), tibble::as_tibble(rc))
  names(row_coord) = c("row", dims)
  col_coord = tibble::tibble(col = rownames(cs), tibble::as_tibble(cs))
  names(col_coord) = c("col", dims)
  inertia_vals = fit$cor^2
  inertia = tibble::tibble(
    dimension = seq_along(inertia_vals),
    inertia = inertia_vals,
    prop = inertia_vals / sum(inertia_vals),
    cum = cumsum(inertia_vals / sum(inertia_vals))
  )
  # Full-table chi-square statistic (independence test on the complete
  # contingency table), independent of the requested number of dimensions:
  # chi2 = N * sum of ALL squared singular correlations of the standardized
  # residual matrix, which equals stats::chisq.test(tab)$statistic.
  E = outer(rowSums(tab), colSums(tab)) / sum(tab)
  sv = svd((tab - E) / sqrt(E))
  chi2_full = sum(sv$d[sv$d > 1e-12]^2)
  meta = tibble::tibble(
    n_used = sum(keep),
    n_dropped = sum(!keep),
    n_row = nr,
    n_col = nc,
    nf = nf,
    chi2_stat = chi2_full
  )
  structure(
    list(
      row_coord = row_coord,
      col_coord = col_coord,
      inertia = inertia,
      meta = meta
    ),
    class = "mv_corresp"
  )
}

#' @export
print.mv_corresp = function(x, ...) {
  cat("Correspondence analysis (mv_corresp)\n")
  cat(sprintf(
    "Table: %d x %d; chi-square statistic: %.3f\n\n",
    x$meta$n_row, x$meta$n_col, x$meta$chi2_stat
  ))
  print(x$inertia)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Multidimensional scaling
# -----------------------------------------------------------------------------

#' Multidimensional Scaling
#'
#' Computes a multidimensional scaling configuration from the distances of the
#' selected numeric columns: classical (metric) MDS via `stats::cmdscale` or
#' non-metric MDS via `MASS::isoMDS`.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to scale.
#' @param k Number of dimensions.
#' @param metric Logical. If `TRUE` (default), classical metric MDS
#'   (`cmdscale`); if `FALSE`, non-metric MDS (`MASS::isoMDS`).
#' @param dist_method Distance method passed to `stats::dist`.
#' @param scale Logical. If `TRUE` (default), columns are standardized first.
#' @param maxit Maximum iterations for `MASS::isoMDS`.
#'
#' @return An object of class `mv_mds`: a list with components `points`
#'   (tibble with `.row` and `Dim1`...), `eigen` (tibble with `dimension` and
#'   `eigenvalue`, metric MDS only), `stress` (single-row tibble, non-metric
#'   MDS only), `gof` (single-row tibble with `gof1`, `gof2`, metric MDS
#'   only), and `meta`.
#'
#' @examples
#' r = mv_mds(mtcars, cyl:carb, k = 2)
#' r$gof
#' plot_mv_mds(r)
#' @export
mv_mds = function(data, cols, k = 2, metric = TRUE,
                  dist_method = "euclidean", scale = TRUE, maxit = 50) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  dist_method = match.arg(dist_method, c(
    "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  ))
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  if (k >= prep$n_used) {
    stop("k must be smaller than the number of complete observations.",
         call. = FALSE)
  }
  if (scale) .mv_check_constant(prep$mat)
  mat = if (scale) .mv_scale(prep$mat) else prep$mat
  d = stats::dist(mat, method = dist_method)
  eigen_tb = NULL
  gof = NULL
  stress_tb = NULL
  if (metric) {
    fit = stats::cmdscale(d, k = k, eig = TRUE)
    pts = fit$points
    colnames(pts) = paste0("Dim", seq_len(k))
    eigen_tb = tibble::tibble(
      dimension = seq_along(fit$eig),
      eigenvalue = fit$eig
    )
    gof = tibble::tibble(gof1 = fit$GOF[1], gof2 = fit$GOF[2])
  } else {
    fit = MASS::isoMDS(d, k = k, maxit = maxit, trace = FALSE)
    pts = fit$points
    colnames(pts) = paste0("Dim", seq_len(k))
    stress_tb = tibble::tibble(stress = fit$stress)
  }
  points = tibble::tibble(.row = prep$row, tibble::as_tibble(pts))
  meta = tibble::tibble(
    n_used = prep$n_used,
    n_dropped = prep$n_dropped,
    k = k,
    metric = metric,
    dist_method = dist_method,
    scale = scale
  )
  structure(
    list(
      points = points,
      eigen = eigen_tb,
      stress = stress_tb,
      gof = gof,
      meta = meta
    ),
    class = "mv_mds"
  )
}

#' @export
print.mv_mds = function(x, ...) {
  cat("Multidimensional scaling (mv_mds)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); k = %d; metric: %s\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$k, x$meta$metric
  ))
  if (!is.null(x$stress)) {
    cat(sprintf("Stress: %.3f\n", x$stress$stress))
  } else {
    cat(sprintf("GOF: %.3f (best), %.3f (second)\n", x$gof$gof1, x$gof$gof2))
  }
  print(utils::head(x$points))
  invisible(x)
}

# -----------------------------------------------------------------------------
# Canonical correlation analysis
# -----------------------------------------------------------------------------

#' Canonical Correlation Analysis
#'
#' Relates two blocks of numeric columns via `stats::cancor` on standardized
#' data, reporting canonical correlations with Wilks' lambda significance
#' tests, raw coefficients, and structure correlations.
#'
#' @param data A data frame.
#' @param x_cols <[`tidy-select`][tidyselect::language]> Numeric columns of the first set.
#' @param y_cols <[`tidy-select`][tidyselect::language]> Numeric columns of the second set;
#'   must not overlap `x_cols`.
#' @param scale Logical. If `TRUE` (default), both sets are standardized so
#'   coefficients live on the standardized scale.
#'
#' @return An object of class `mv_cancor`: a list with components
#'   `correlations` (tibble with `dimension`, `canonical_correlation`,
#'   `wilks_lambda`, `chi_square`, `df`, `p_value`), `x_coef` (tibble with
#'   `variable` and `Can1`...), `y_coef`, `x_structure` (tibble with
#'   `variable` and `Can1`... structure correlations), `y_structure`, and
#'   `meta` (where `n_dropped` counts rows incomplete in *either* column
#'   set).
#'
#' @examples
#' r = mv_cancor(mtcars, mpg:disp, hp:wt)
#' r$correlations
#' @export
mv_cancor = function(data, x_cols, y_cols, scale = TRUE) {
  .mv_frame(data)
  px = .mv_matrix(dplyr::select(data, {{ x_cols }}), arg = "x_cols")
  py = .mv_matrix(dplyr::select(data, {{ y_cols }}), arg = "y_cols")
  overlap = intersect(px$names, py$names)
  if (length(overlap) > 0) {
    stop(
      "x_cols and y_cols must not share column(s): ",
      paste(overlap, collapse = ", "), ".",
      call. = FALSE
    )
  }
  rows = intersect(px$row, py$row)
  if (length(rows) < 3) {
    stop("Need at least 3 complete observations across both column sets.",
         call. = FALSE)
  }
  X = px$mat[match(rows, px$row), , drop = FALSE]
  Y = py$mat[match(rows, py$row), , drop = FALSE]
  p = ncol(X)
  q = ncol(Y)
  n = length(rows)
  if (scale) {
    .mv_check_constant(X)
    .mv_check_constant(Y)
  }
  Xs = if (scale) .mv_scale(X) else X
  Ys = if (scale) .mv_scale(Y) else Y
  if (n <= p + q) {
    stop("Need more complete observations than p + q selected columns.",
         call. = FALSE)
  }
  cc = stats::cancor(Xs, Ys)
  r = length(cc$cor)
  dims = paste0("Can", seq_len(r))
  wilks = chi2 = dfs = pval = numeric(r)
  for (i in seq_len(r)) {
    lam = prod(1 - cc$cor[i:r]^2)
    wilks[i] = lam
    chi2[i] = -(n - 1 - (p + q + 1) / 2) * log(lam)
    dfs[i] = (p - i + 1) * (q - i + 1)
    pval[i] = stats::pchisq(chi2[i], df = dfs[i], lower.tail = FALSE)
  }
  correlations = tibble::tibble(
    dimension = seq_len(r),
    canonical_correlation = cc$cor,
    wilks_lambda = wilks,
    chi_square = chi2,
    df = dfs,
    p_value = pval
  )
  xc = cc$xcoef
  colnames(xc) = dims
  yc = cc$ycoef
  colnames(yc) = dims
  x_coef = tibble::tibble(variable = colnames(Xs), tibble::as_tibble(xc))
  names(x_coef) = c("variable", dims)
  y_coef = tibble::tibble(variable = colnames(Ys), tibble::as_tibble(yc))
  names(y_coef) = c("variable", dims)
  U = Xs %*% cc$xcoef
  V = Ys %*% cc$ycoef
  colnames(U) = dims
  colnames(V) = dims
  xs_str = stats::cor(Xs, U)
  colnames(xs_str) = dims
  ys_str = stats::cor(Ys, V)
  colnames(ys_str) = dims
  x_structure = tibble::tibble(variable = colnames(Xs), tibble::as_tibble(xs_str))
  names(x_structure) = c("variable", dims)
  y_structure = tibble::tibble(variable = colnames(Ys), tibble::as_tibble(ys_str))
  names(y_structure) = c("variable", dims)
  meta = tibble::tibble(
    n_used = n,
    n_dropped = nrow(data) - n,
    p = p,
    q = q,
    r = r,
    scale = scale
  )
  structure(
    list(
      correlations = correlations,
      x_coef = x_coef,
      y_coef = y_coef,
      x_structure = x_structure,
      y_structure = y_structure,
      meta = meta
    ),
    class = "mv_cancor"
  )
}

#' @export
print.mv_cancor = function(x, ...) {
  cat("Canonical correlation analysis (mv_cancor)\n")
  cat(sprintf(
    "Observations: %d used; x variables: %d; y variables: %d\n\n",
    x$meta$n_used, x$meta$p, x$meta$q
  ))
  print(x$correlations)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Partial correlation
# -----------------------------------------------------------------------------

#' Partial Correlation Matrix
#'
#' Computes first-order-controlled partial correlations among the selected
#' numeric columns by inverting the correlation matrix, with t-tests for each
#' pairwise partial correlation. Each partial correlation controls for the
#' other `p - 2` variables and is tested with `df = n - p` degrees of
#' freedom.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns; the set must
#'   not be (near-)singular.
#'
#' @return An object of class `mv_pcor`: a list with components
#'   `partial_corr` (long tibble with `variable_1`, `variable_2`,
#'   `partial_corr`, `t_stat`, `df`, `p_value`), `pcor_matrix` (wide tibble
#'   with `variable` and one column per variable), and `meta`.
#'
#' @examples
#' r = mv_pcor(mtcars, mpg:wt)
#' r$partial_corr
#' @export
mv_pcor = function(data, cols) {
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  n = prep$n_used
  p = ncol(prep$mat)
  if (p < 3) {
    stop("Partial correlation needs at least 3 variables (2 to control).",
         call. = FALSE)
  }
  if (n <= p + 2) {
    stop("Need more complete observations than p + 2 for partial correlation tests.",
         call. = FALSE)
  }
  S = stats::cor(prep$mat)
  Rinv = tryCatch(
    solve(S),
    error = function(e) {
      stop(
        "The correlation matrix is singular; drop redundant (collinear) variables.",
        call. = FALSE
      )
    }
  )
  dvec = sqrt(diag(Rinv))
  P = -Rinv / outer(dvec, dvec)
  diag(P) = 1
  prs = utils::combn(p, 2)
  df = n - p
  r_vec = numeric(ncol(prs))
  t_vec = numeric(ncol(prs))
  p_vec = numeric(ncol(prs))
  for (j in seq_len(ncol(prs))) {
    i1 = prs[1, j]
    i2 = prs[2, j]
    rj = P[i1, i2]
    r_vec[j] = rj
    t_vec[j] = rj * sqrt(df / (1 - rj^2))
    p_vec[j] = 2 * stats::pt(-abs(t_vec[j]), df = df)
  }
  partial_corr = tibble::tibble(
    variable_1 = prep$names[prs[1, ]],
    variable_2 = prep$names[prs[2, ]],
    partial_corr = r_vec,
    t_stat = t_vec,
    df = df,
    p_value = p_vec
  )
  colnames(P) = prep$names
  rownames(P) = prep$names
  pcor_matrix = tibble::tibble(variable = prep$names, tibble::as_tibble(P))
  meta = tibble::tibble(
    n_used = n,
    n_dropped = prep$n_dropped,
    p = p,
    df = df
  )
  structure(
    list(
      partial_corr = partial_corr,
      pcor_matrix = pcor_matrix,
      meta = meta
    ),
    class = "mv_pcor"
  )
}

#' @export
print.mv_pcor = function(x, ...) {
  cat("Partial correlations (mv_pcor)\n")
  cat(sprintf(
    "Observations: %d used, %d dropped (NA); variables: %d; df: %d\n\n",
    x$meta$n_used, x$meta$n_dropped, x$meta$p, x$meta$df
  ))
  print(x$partial_corr)
  invisible(x)
}

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

#' Scree Plot
#'
#' Plots the proportion and cumulative proportion of variance explained by
#' each component (PCA) or factor (EFA).
#'
#' @param x An `mv_pca` or `mv_efa` result.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' r = mv_pca(mtcars, cyl:carb)
#' plot_mv_scree(r)
#' @export
plot_mv_scree = function(x) {
  if (!inherits(x, c("mv_pca", "mv_efa"))) {
    stop("x must be an mv_pca or mv_efa result.", call. = FALSE)
  }
  v = x$variance
  comp_name = names(v)[1]
  long = tibble::tibble(
    component = rep(factor(v[[1]], levels = v[[1]]), 2),
    measure = rep(c("Proportion", "Cumulative"), each = nrow(v)),
    value = c(v$prop_var, v$cum_var)
  )
  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$component, y = .data$value,
    group = .data$measure, color = .data$measure
  )) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = comp_name,
      y = "Proportion of variance",
      color = NULL
    ) +
    ggplot2::theme_minimal()
}

#' Scores Plot
#'
#' Plots the first two score dimensions of an `mv_pca` result (PC1 vs PC2) or
#' an LDA `mv_discrim` result (LD1 vs LD2), optionally colored by group.
#'
#' @param x An `mv_pca` result, or an `mv_discrim` result from LDA (QDA has
#'   no scores); LDA with only 2 groups yields a single discriminant
#'   function, which is not enough for a 2-D scores plot.
#' @param group Optional vector with one value per row of the scores, used
#'   for coloring.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' r = mv_pca(mtcars, cyl:carb)
#' plot_mv_scores(r)
#' @export
plot_mv_scores = function(x, group = NULL) {
  if (inherits(x, "mv_pca")) {
    sc = x$scores
    if (ncol(sc) < 3) {
      stop("PCA needs at least 2 components for a scores plot.", call. = FALSE)
    }
    c1 = "PC1"
    c2 = "PC2"
    l1 = sprintf("PC1 (%.1f%%)", 100 * x$variance$prop_var[1])
    l2 = sprintf("PC2 (%.1f%%)", 100 * x$variance$prop_var[2])
  } else if (inherits(x, "mv_discrim")) {
    if (is.null(x$scores)) {
      stop("QDA results have no scores; use an LDA result.", call. = FALSE)
    }
    sc = x$scores
    if (ncol(sc) < 3) {
      stop(
        "LDA scores plot needs at least 2 discriminant functions, ",
        "i.e. at least 3 groups.",
        call. = FALSE
      )
    }
    c1 = names(sc)[2]
    c2 = names(sc)[3]
    l1 = c1
    l2 = c2
  } else {
    stop("x must be an mv_pca or mv_discrim (LDA) result.", call. = FALSE)
  }
  if (!is.null(group)) {
    if (length(group) != nrow(sc)) {
      stop("group must have one value per row of the scores.", call. = FALSE)
    }
    sc$group = factor(group)
  }
  p = ggplot2::ggplot(
    sc,
    ggplot2::aes(x = .data[[c1]], y = .data[[c2]])
  )
  if (!is.null(group)) {
    p = p + ggplot2::geom_point(ggplot2::aes(color = .data$group))
  } else {
    p = p + ggplot2::geom_point(color = "steelblue")
  }
  p +
    ggplot2::labs(x = l1, y = l2, color = NULL) +
    ggplot2::theme_minimal()
}

#' Dendrogram Plot
#'
#' Draws the dendrogram of an `mv_hclust` result with ggplot2, with a dashed
#' reference line at the tree height where the `k` clusters are cut.
#'
#' @param x An `mv_hclust` result.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' r = mv_hclust(mtcars, cyl:carb, k = 3)
#' plot_mv_dendrogram(r)
#' @export
plot_mv_dendrogram = function(x) {
  if (!inherits(x, "mv_hclust")) {
    stop("x must be an mv_hclust result.", call. = FALSE)
  }
  dd = .mv_dendro_segments(x$hclust)
  # Cut height producing meta$k clusters: the height of the (n - k)-th merge.
  n = length(x$hclust$order)
  cut_h = x$hclust$height[n - x$meta$k]
  ggplot2::ggplot(dd$segments, ggplot2::aes(
    x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend
  )) +
    ggplot2::geom_segment() +
    ggplot2::geom_hline(
      yintercept = cut_h, linetype = "dashed", color = "grey40"
    ) +
    ggplot2::scale_x_continuous(
      breaks = dd$labels$x, labels = dd$labels$label
    ) +
    ggplot2::labs(
      x = NULL, y = "Height",
      subtitle = sprintf("Cut height for k = %d: %.3f", x$meta$k, cut_h)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}

#' Silhouette Plot
#'
#' Draws the silhouette plot (per-observation silhouette widths grouped by
#' cluster, sorted within clusters) with the average width as a reference
#' line. Accepts an `mv_cluster_metrics`, `mv_kmeans`, or `mv_hclust` result;
#' for the latter two the widths are computed from the stored analysis matrix,
#' using the distance method recorded in the result's `meta` (euclidean for
#' `mv_kmeans`, which always clusters in euclidean space).
#'
#' @param x An `mv_cluster_metrics`, `mv_kmeans`, or `mv_hclust` result.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' km = mv_kmeans(iris, Sepal.Length:Petal.Width, centers = 3)
#' plot_mv_silhouette(km)
#' @export
plot_mv_silhouette = function(x) {
  if (inherits(x, "mv_cluster_metrics")) {
    sil = x$silhouette
  } else if (inherits(x, "mv_kmeans")) {
    sil = .mv_sil_df(x$data_mat, x$cluster$cluster)
  } else if (inherits(x, "mv_hclust")) {
    sil = .mv_sil_df(
      x$data_mat, x$cluster$cluster,
      dist_method = x$meta$dist_method
    )
  } else {
    stop(
      "x must be an mv_cluster_metrics, mv_kmeans, or mv_hclust result.",
      call. = FALSE
    )
  }
  ord = order(as.integer(sil$cluster), -sil$sil_width)
  sil = sil[ord, , drop = FALSE]
  sil$pos = seq_len(nrow(sil))
  avg = mean(sil$sil_width)
  ggplot2::ggplot(sil, ggplot2::aes(
    x = .data$pos, y = .data$sil_width, fill = .data$cluster
  )) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(
      yintercept = avg, linetype = "dashed", color = "grey30"
    ) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::labs(
      x = "Observation (ordered within cluster)",
      y = "Silhouette width",
      fill = "Cluster",
      subtitle = sprintf("Average silhouette width = %.3f", avg)
    ) +
    ggplot2::theme_minimal()
}

#' Elbow Plot
#'
#' Plots the total within-cluster sum of squares against the number of
#' clusters `k`, for k-means (refit per `k`) or hierarchical clustering
#' (one tree, cut at each `k`). The k-means branch fixes the RNG seed
#' (1234) internally so the curve is reproducible; the caller's RNG state
#' is restored afterwards.
#'
#' @param data A data frame.
#' @param cols <[`tidy-select`][tidyselect::language]> Numeric columns to cluster on.
#' @param max_k Largest number of clusters to evaluate.
#' @param method Either `"kmeans"` (default) or `"hclust"`.
#' @param nstart Number of random starts for k-means.
#' @param dist_method Distance method for the hierarchical tree.
#' @param hclust_method Agglomeration method for the hierarchical tree.
#' @param scale Logical. If `TRUE` (default), columns are standardized first.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' plot_mv_elbow(mtcars, cyl:carb, max_k = 6)
#' plot_mv_elbow(mtcars, cyl:carb, max_k = 6, method = "hclust")
#' @export
plot_mv_elbow = function(data, cols, max_k = 10,
                         method = c("kmeans", "hclust"), nstart = 25,
                         dist_method = "euclidean",
                         hclust_method = "ward.D2", scale = TRUE) {
  method = match.arg(method)
  .mv_frame(data)
  prep = .mv_matrix(dplyr::select(data, {{ cols }}))
  if (scale) .mv_check_constant(prep$mat)
  mat = if (scale) .mv_scale(prep$mat) else prep$mat
  n = nrow(mat)
  if (!is.numeric(max_k) || length(max_k) != 1 || max_k < 2 ||
      max_k != floor(max_k)) {
    stop("max_k must be an integer >= 2.", call. = FALSE)
  }
  max_k = min(max_k, n - 1)
  ks = seq_len(max_k)
  wss = if (method == "kmeans") {
    # Fix the RNG only for this branch and restore the caller's seed state.
    had_seed = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    seed0 = if (had_seed) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    on.exit(
      if (had_seed) {
        assign(".Random.seed", seed0, envir = globalenv(), inherits = FALSE)
      },
      add = TRUE
    )
    set.seed(1234)
    vapply(ks, function(k) {
      stats::kmeans(mat, centers = k, nstart = nstart)$tot.withinss
    }, numeric(1))
  } else {
    hc = stats::hclust(
      stats::dist(mat, method = dist_method),
      method = hclust_method
    )
    vapply(ks, function(k) {
      .mv_wss(mat, stats::cutree(hc, k = k))
    }, numeric(1))
  }
  df = tibble::tibble(k = ks, wss = wss)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$k, y = .data$wss)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_point(color = "steelblue") +
    ggplot2::labs(
      x = "Number of clusters (k)",
      y = "Total within-cluster sum of squares"
    ) +
    ggplot2::theme_minimal()
}

#' Correspondence Analysis Map
#'
#' Draws the symmetric two-dimensional map of a `mv_corresp` result: row and
#' column coordinates with labels.
#'
#' @param x An `mv_corresp` result with at least two dimensions.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(42)
#' df = data.frame(
#'   breed = rep(c("A", "B", "C"), each = 30),
#'   grade = sample(c("low", "mid", "high"), 90, TRUE)
#' )
#' r = mv_corresp(df, breed, grade)
#' plot_mv_corresp(r)
#' @export
plot_mv_corresp = function(x) {
  if (!inherits(x, "mv_corresp")) {
    stop("x must be an mv_corresp result.", call. = FALSE)
  }
  if (ncol(x$row_coord) < 3) {
    stop("Correspondence analysis needs nf >= 2 for a 2-D map.",
         call. = FALSE)
  }
  xs = c(x$row_coord$Dim1, x$col_coord$Dim1)
  ys = c(x$row_coord$Dim2, x$col_coord$Dim2)
  labs = c(x$row_coord$row, x$col_coord$col)
  grp = rep(c("Row", "Column"), each = nrow(x$row_coord))
  df = tibble::tibble(
    x = xs,
    y = ys,
    label = labs,
    set = grp
  )
  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$x, y = .data$y, color = .data$set, shape = .data$set
  )) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      vjust = -0.8, size = 3, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    ggplot2::labs(x = "Dim1", y = "Dim2", color = NULL, shape = NULL) +
    ggplot2::theme_minimal()
}

#' MDS Configuration Plot
#'
#' Plots the first two dimensions of an `mv_mds` configuration.
#'
#' @param x An `mv_mds` result with `k >= 2`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' r = mv_mds(mtcars, cyl:carb, k = 2)
#' plot_mv_mds(r)
#' @export
plot_mv_mds = function(x) {
  if (!inherits(x, "mv_mds")) {
    stop("x must be an mv_mds result.", call. = FALSE)
  }
  if (!all(c("Dim1", "Dim2") %in% names(x$points))) {
    stop("MDS needs k >= 2 for a 2-D configuration.", call. = FALSE)
  }
  sub = if (!is.null(x$stress)) {
    sprintf("Stress = %.3f", x$stress$stress)
  } else {
    sprintf("GOF = %.3f", x$gof$gof1)
  }
  ggplot2::ggplot(
    x$points,
    ggplot2::aes(x = .data$Dim1, y = .data$Dim2)
  ) +
    ggplot2::geom_point(color = "steelblue") +
    ggplot2::labs(x = "Dim1", y = "Dim2", subtitle = sub) +
    ggplot2::theme_minimal()
}

#' @importFrom rlang .data
NULL
