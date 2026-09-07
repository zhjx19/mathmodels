# Tests for the multivariate statistics toolkit (mv_stats.R)

mt_cols = c("cyl", "disp", "hp", "drat", "wt", "qsec", "vs", "am", "gear", "carb")
mt = mtcars[, mt_cols]

# ---------------------------------------------------------------------------
# mv_pca
# ---------------------------------------------------------------------------

test_that("mv_pca matches prcomp", {
  r = mv_pca(mtcars, cyl:carb)
  ref = stats::prcomp(mt, center = TRUE, scale. = TRUE)
  expect_s3_class(r, "mv_pca")
  expect_s3_class(r$scores, "tbl_df")
  expect_s3_class(r$loadings, "tbl_df")
  expect_equal(r$variance$sdev, ref$sdev, tolerance = 1e-10)
  expect_equal(
    as.matrix(abs(r$scores[, -1])),
    abs(ref$x),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_equal(
    as.matrix(abs(r$loadings[, -1])),
    abs(unclass(ref$rotation)),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_equal(sum(r$variance$prop_var), 1, tolerance = 1e-10)
  expect_equal(
    r$variance$cum_var,
    cumsum(r$variance$prop_var),
    tolerance = 1e-10
  )
})

test_that("mv_pca drops NA rows and reports counts", {
  df = mtcars[, c("mpg", "disp", "hp")]
  df$mpg[3] = NA
  r = mv_pca(df, mpg:hp)
  expect_equal(r$meta$n_used, 31)
  expect_equal(r$meta$n_dropped, 1)
  expect_equal(nrow(r$scores), 31)
})

test_that("mv_pca input validation", {
  df = data.frame(a = 1:4, b = c("x", "y", "z", "w"))
  expect_error(mv_pca(df, a:b), "numeric")
  expect_error(mv_pca(mtcars, nonexistent))
})

test_that("mv_pca rejects constant columns when scaling", {
  df = mtcars[, c("mpg", "disp")]
  df$const = 5
  expect_error(mv_pca(df, mpg:const), "constant")
})

# ---------------------------------------------------------------------------
# mv_efa
# ---------------------------------------------------------------------------

test_that("mv_efa returns tidy components", {
  r = mv_efa(mtcars, cyl:carb, factors = 2)
  expect_s3_class(r, "mv_efa")
  expect_s3_class(r$loadings, "tbl_df")
  expect_equal(names(r$loadings), c("variable", "Factor1", "Factor2", "uniqueness"))
  expect_true(all(r$loadings$uniqueness >= 0 & r$loadings$uniqueness <= 1))
  expect_equal(nrow(r$scores), r$meta$n_used)
  expect_equal(ncol(r$scores), 3) # .row + 2 factors
  expect_equal(sum(r$variance$prop_var), sum(r$variance$ss_loadings) / 10)
  expect_equal(
    r$variance$cum_var,
    cumsum(r$variance$prop_var),
    tolerance = 1e-10
  )
})

test_that("mv_efa can skip scores", {
  r = mv_efa(mtcars, cyl:carb, factors = 2, scores = "none")
  expect_null(r$scores)
})

test_that("mv_efa input validation", {
  expect_error(mv_efa(mtcars, cyl:carb, factors = 10), "smaller")
  expect_error(mv_efa(mtcars, cyl:carb, factors = 1.5), "integer")
})

# ---------------------------------------------------------------------------
# mv_hclust
# ---------------------------------------------------------------------------

test_that("mv_hclust matches hclust + cutree", {
  r = mv_hclust(mtcars, cyl:carb, k = 3)
  d = stats::dist(scale(mt))
  h = stats::hclust(d, method = "ward.D2")
  expect_s3_class(r, "mv_hclust")
  expect_equal(
    as.integer(r$cluster$cluster),
    unname(stats::cutree(h, k = 3))
  )
  expect_s3_class(r$hclust, "hclust")
  expect_equal(max(r$heights$height), max(h$height), tolerance = 1e-10)
  expect_equal(nrow(r$heights), 31)
  expect_equal(
    unname(apply(r$data_mat, 2, stats::sd)),
    rep(1, 10),
    tolerance = 1e-10
  )
})

test_that("mv_hclust without scaling matches raw-data clustering", {
  r = mv_hclust(mtcars, cyl:carb, k = 2, scale = FALSE)
  h = stats::hclust(stats::dist(mt), method = "ward.D2")
  expect_equal(
    as.integer(r$cluster$cluster),
    unname(stats::cutree(h, k = 2))
  )
})

test_that("mv_hclust input validation", {
  expect_error(mv_hclust(mtcars, cyl:carb, k = 32), "smaller")
  expect_error(mv_hclust(mtcars, cyl:carb, method = "nope"), "should be one of")
})

# ---------------------------------------------------------------------------
# mv_kmeans
# ---------------------------------------------------------------------------

test_that("mv_kmeans matches kmeans", {
  set.seed(1)
  r = mv_kmeans(mtcars, cyl:carb, centers = 3)
  set.seed(1)
  ref = stats::kmeans(scale(mt), centers = 3, nstart = 25)
  expect_s3_class(r, "mv_kmeans")
  expect_equal(as.integer(r$cluster$cluster), as.integer(ref$cluster))
  expect_equal(r$meta$tot_withinss, ref$tot.withinss, tolerance = 1e-10)
  expect_equal(sum(r$within_ss$size), 32)
  expect_equal(sum(r$within_ss$share), 1, tolerance = 1e-10)
})

test_that("mv_kmeans input validation", {
  expect_error(mv_kmeans(mtcars, cyl:carb, centers = 3, nstart = 0), "nstart")
})

# ---------------------------------------------------------------------------
# mv_cluster_metrics
# ---------------------------------------------------------------------------

test_that("silhouette widths match cluster::silhouette", {
  skip_if_not_installed("cluster")
  set.seed(7)
  df = data.frame(
    a = c(rnorm(20), rnorm(20, mean = 5)),
    b = c(rnorm(20), rnorm(20, mean = 5))
  )
  cl = rep(c("g1", "g2"), each = 20)
  r = mv_cluster_metrics(df, a:b, cl, scale = FALSE)
  ref = cluster::silhouette(rep(1:2, each = 20), stats::dist(df))
  expect_s3_class(r, "mv_cluster_metrics")
  expect_equal(
    as.numeric(r$silhouette$sil_width),
    as.numeric(ref[, "sil_width"]),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(r$summary$value[r$summary$metric == "avg_silhouette"]),
    as.numeric(mean(ref[, "sil_width"])),
    tolerance = 1e-8
  )
})

test_that("validity metrics rank separated clusters above random ones", {
  set.seed(11)
  good = data.frame(
    a = c(rnorm(30), rnorm(30, mean = 8)),
    b = c(rnorm(30), rnorm(30, mean = 8))
  )
  labels = rep(c("x", "y"), each = 30)
  r_good = mv_cluster_metrics(good, a:b, labels)
  r_rand = mv_cluster_metrics(good, a:b, sample(labels))
  ch = r_good$summary$value[r_good$summary$metric == "calinski_harabasz"]
  ch_r = r_rand$summary$value[r_rand$summary$metric == "calinski_harabasz"]
  db = r_good$summary$value[r_good$summary$metric == "davies_bouldin"]
  db_r = r_rand$summary$value[r_rand$summary$metric == "davies_bouldin"]
  expect_gt(ch, ch_r)
  expect_lt(db, db_r)
  expect_gt(
    r_good$summary$value[r_good$summary$metric == "dunn"], 1
  )
  expect_true(all(r_good$silhouette$sil_width <= 1))
  expect_true(all(r_good$silhouette$sil_width >= -1))
})

test_that("mv_cluster_metrics input validation", {
  expect_error(
    mv_cluster_metrics(mtcars, cyl:carb, rep(1, 32)),
    "at least 2"
  )
  expect_error(
    mv_cluster_metrics(mtcars, cyl:carb, rep(1:2, 16)[-1]),
    "one label per row"
  )
  expect_error(
    mv_cluster_metrics(mtcars, cyl:carb, c(NA, rep(1:2, each = 15), 1)),
    "NA"
  )
})

# ---------------------------------------------------------------------------
# mv_discrim
# ---------------------------------------------------------------------------

test_that("mv_discrim lda matches MASS::lda", {
  r = mv_discrim(iris, Species, Sepal.Length:Petal.Width)
  ref = MASS::lda(x = iris[, 1:4], grouping = iris$Species)
  pred = stats::predict(ref, newdata = iris[, 1:4])
  expect_s3_class(r, "mv_discrim")
  expect_equal(
    as.character(r$predictions$predicted),
    as.character(pred$class)
  )
  expect_equal(r$meta$accuracy, mean(pred$class == iris$Species))
  expect_equal(sum(r$confusion$n), 150)
  expect_s3_class(r$scores, "tbl_df")
  expect_equal(ncol(r$scores), 3) # .row + LD1 + LD2
  expect_equal(names(r$structure), c("variable", "LD1", "LD2"))
})

test_that("mv_discrim qda has no scores and reports means", {
  r = mv_discrim(iris, Species, Sepal.Length:Petal.Width, method = "qda")
  expect_null(r$scores)
  expect_equal(nrow(r$structure), 3)
  expect_s3_class(r$structure$group, "factor")
})

test_that("mv_discrim classifies newdata", {
  r = mv_discrim(iris, Species, Sepal.Length:Petal.Width)
  nd = iris[1:5, 1:4]
  np = r$new_predictions
  expect_null(np)
  r2 = mv_discrim(iris, Species, Sepal.Length:Petal.Width, newdata = nd)
  expect_equal(nrow(r2$new_predictions), 5)
  expect_true(all(as.character(r2$new_predictions$predicted) %in%
    levels(iris$Species)))
})

test_that("mv_discrim input validation", {
  expect_error(
    mv_discrim(iris, Species, Species:Petal.Width),
    "must not be among the predictors"
  )
  df = iris
  df$Species[1] = NA
  expect_error(mv_discrim(df, Species, Sepal.Length:Petal.Width), "NA")
  df2 = data.frame(g = rep("only", 10), x = 1:10, y = rnorm(10))
  expect_error(mv_discrim(df2, g, x:y), "at least 2 levels")
})

# ---------------------------------------------------------------------------
# mv_corresp
# ---------------------------------------------------------------------------

test_that("mv_corresp matches MASS::corresp", {
  set.seed(42)
  df = data.frame(
    breed = rep(c("A", "B", "C"), each = 30),
    grade = sample(c("low", "mid", "high"), 90, TRUE)
  )
  r = mv_corresp(df, breed, grade, nf = 2)
  f = MASS::corresp(unclass(table(df$breed, df$grade)), nf = 2)
  expect_s3_class(r, "mv_corresp")
  expect_equal(r$row_coord$Dim1, as.numeric(f$rscore[, 1]), tolerance = 1e-8)
  expect_equal(r$col_coord$Dim2, as.numeric(f$cscore[, 2]), tolerance = 1e-8)
  expect_equal(
    sum(r$inertia$inertia),
    sum(f$cor^2),
    tolerance = 1e-8
  )
  expect_equal(
    r$meta$chi2_stat,
    sum(table(df$breed, df$grade)) * sum(f$cor^2),
    tolerance = 1e-8
  )
  expect_equal(sum(r$inertia$prop), 1, tolerance = 1e-10)
})

test_that("mv_corresp input validation", {
  set.seed(42)
  df = data.frame(
    breed = rep(c("A", "B", "C"), each = 30),
    grade = sample(c("low", "mid", "high"), 90, TRUE)
  )
  expect_error(mv_corresp(df, breed, breed), "different")
  expect_error(mv_corresp(df, breed, grade, nf = 5), "smaller")
  expect_error(mv_corresp(df, breed, grade, nf = 0), "positive")
})

# ---------------------------------------------------------------------------
# mv_mds
# ---------------------------------------------------------------------------

test_that("mv_mds metric matches cmdscale up to reflection", {
  r = mv_mds(mtcars, cyl:carb, k = 2, scale = FALSE)
  f = stats::cmdscale(stats::dist(mt), k = 2, eig = TRUE)
  expect_s3_class(r, "mv_mds")
  expect_equal(
    stats::dist(as.matrix(r$points[, c("Dim1", "Dim2")])),
    stats::dist(f$points),
    ignore_attr = TRUE
  )
  expect_true(all(r$gof$gof1 >= 0 & r$gof$gof1 <= 1))
  expect_equal(ncol(r$eigen), 2)
})

test_that("mv_mds non-metric reports stress", {
  r = mv_mds(mtcars, cyl:carb, k = 2, metric = FALSE)
  expect_null(r$eigen)
  expect_null(r$gof)
  expect_equal(nrow(r$stress), 1)
  expect_true(r$stress$stress >= 0 && r$stress$stress <= 100)
})

test_that("mv_mds input validation", {
  expect_error(mv_mds(mtcars, cyl:carb, k = 32), "smaller")
})

# ---------------------------------------------------------------------------
# mv_cancor
# ---------------------------------------------------------------------------

test_that("mv_cancor matches cancor", {
  x_cols = c("mpg", "cyl", "disp")
  y_cols = c("hp", "drat", "wt")
  r = mv_cancor(mtcars, all_of(x_cols), all_of(y_cols))
  cc = stats::cancor(scale(mtcars[, x_cols]), scale(mtcars[, y_cols]))
  expect_s3_class(r, "mv_cancor")
  expect_equal(
    r$correlations$canonical_correlation,
    cc$cor,
    tolerance = 1e-10
  )
  expect_true(all(r$correlations$p_value >= 0 & r$correlations$p_value <= 1))
  expect_equal(r$meta$r, 3)
  expect_equal(nrow(r$x_structure), 3)
  expect_equal(nrow(r$y_structure), 3)
})

test_that("mv_cancor input validation", {
  expect_error(
    mv_cancor(mtcars, mpg:disp, disp:wt),
    "must not share"
  )
  df = data.frame(x1 = 1:3, x2 = 2:4, y1 = 3:5)
  expect_error(mv_cancor(df, x1:x2, y1), "p \\+ q")
})

# ---------------------------------------------------------------------------
# mv_pcor
# ---------------------------------------------------------------------------

test_that("mv_pcor equals correlation of residuals", {
  df = mtcars[, c("mpg", "disp", "wt")]
  r = mv_pcor(df, mpg:wt)
  rx = stats::residuals(stats::lm(mpg ~ disp, data = df))
  rw = stats::residuals(stats::lm(wt ~ disp, data = df))
  expected = stats::cor(rx, rw)
  i = which(r$partial_corr$variable_1 == "mpg" &
    r$partial_corr$variable_2 == "wt")
  expect_s3_class(r, "mv_pcor")
  expect_equal(r$partial_corr$partial_corr[i], expected, tolerance = 1e-8)
  expect_equal(diag(as.matrix(r$pcor_matrix[, -1])), rep(1, 3),
    tolerance = 1e-10
  )
  expect_equal(r$meta$df, 29) # n - p = 32 - 3 (df = n - 2 - q, q = p - 2)
  # t-statistic must match the residual-correlation t-test with df = n - p
  i2 = which(r$partial_corr$variable_1 == "mpg" &
    r$partial_corr$variable_2 == "wt")
  expect_equal(
    r$partial_corr$t_stat[i2],
    r$partial_corr$partial_corr[i2] * sqrt(
      29 / (1 - r$partial_corr$partial_corr[i2]^2)
    ),
    tolerance = 1e-10
  )
})

test_that("mv_pcor input validation", {
  expect_error(mv_pcor(mtcars, c(mpg, disp)), "at least 3")
  df = data.frame(x = 1:20, y = 2 * (1:20), z = (1:20)^2 + rnorm(20))
  expect_error(mv_pcor(df, x:z), "singular")
})

# ---------------------------------------------------------------------------
# Print methods
# ---------------------------------------------------------------------------

test_that("print methods return the object invisibly", {
  objs = list(
    mv_pca(mtcars, cyl:carb),
    mv_efa(mtcars, cyl:carb, factors = 2),
    mv_hclust(mtcars, cyl:carb, k = 3),
    mv_kmeans(mtcars, cyl:carb, centers = 3),
    mv_cluster_metrics(mtcars, cyl:carb, rep(1:2, each = 16)),
    mv_discrim(iris, Species, Sepal.Length:Petal.Width),
    mv_mds(mtcars, cyl:carb, k = 2),
    mv_cancor(mtcars, mpg:disp, hp:wt),
    mv_pcor(mtcars, mpg:wt)
  )
  for (o in objs) {
    expect_invisible(suppressWarnings(print(o)))
  }
})

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

test_that("plot functions return ggplot objects", {
  pca = mv_pca(mtcars, cyl:carb)
  efa = mv_efa(mtcars, cyl:carb, factors = 2)
  hc = mv_hclust(mtcars, cyl:carb, k = 3)
  km = mv_kmeans(mtcars, cyl:carb, centers = 3)
  met = mv_cluster_metrics(mtcars, cyl:carb, rep(1:2, each = 16))
  lda = mv_discrim(iris, Species, Sepal.Length:Petal.Width)
  cor_ = mv_corresp(
    data.frame(
      g = rep(c("a", "b", "c"), each = 20),
      h = rep(c("x", "y", "z"), 20)
    ),
    g, h
  )
  mds = mv_mds(mtcars, cyl:carb, k = 2)
  plots = list(
    plot_mv_scree(pca),
    plot_mv_scree(efa),
    plot_mv_scores(pca),
    plot_mv_scores(pca, group = rep(c("A", "B"), 16)),
    plot_mv_scores(lda, group = iris$Species),
    plot_mv_dendrogram(hc),
    plot_mv_silhouette(km),
    plot_mv_silhouette(hc),
    plot_mv_silhouette(met),
    plot_mv_elbow(mtcars, cyl:carb, max_k = 4),
    plot_mv_elbow(mtcars, cyl:carb, max_k = 4, method = "hclust"),
    plot_mv_corresp(cor_),
    plot_mv_mds(mds)
  )
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    layer_rows = vapply(ggplot2::ggplot_build(p)$data, nrow, numeric(1))
    expect_true(any(layer_rows > 0))
  }
})

test_that("plot functions validate inputs", {
  expect_error(plot_mv_scree(mtcars), "mv_pca or mv_efa")
  expect_error(plot_mv_dendrogram(mtcars), "mv_hclust")
  expect_error(plot_mv_silhouette(mtcars), "mv_cluster_metrics")
  expect_error(plot_mv_corresp(mtcars), "mv_corresp")
  expect_error(plot_mv_mds(mtcars), "mv_mds")
  expect_error(
    plot_mv_scores(mv_pca(mtcars, cyl:carb), group = 1:3),
    "one value per row"
  )
})

# ---------------------------------------------------------------------------
# Regression tests for review fixes
# ---------------------------------------------------------------------------

test_that("plot_mv_scores rejects 2-group LDA (single LD) with a clear error", {
  ir2 = iris[iris$Species != "setosa", ]
  lda2 = mv_discrim(ir2, Species, Sepal.Length:Petal.Width)
  expect_error(plot_mv_scores(lda2), "at least 3 groups")
  # and no lingering NA aesthetics even before the guard triggers
  expect_error(plot_mv_scores(lda2, group = ir2$Species), "at least 3 groups")
})

test_that("plot_mv_dendrogram keeps labels for tibble input", {
  tb = tibble::as_tibble(mtcars)
  r = mv_hclust(tb, cyl:carb, k = 3)
  labs = .mv_dendro_segments(r$hclust)$labels$label
  expect_false(any(is.na(labs)))
  expect_equal(length(labs), 32)
  # a tibble built from scratch (no rownames at all)
  tb2 = tibble::tibble(a = rnorm(12), b = rnorm(12), c = rnorm(12))
  r2 = mv_hclust(tb2, a:c, k = 2)
  labs2 = .mv_dendro_segments(r2$hclust)$labels$label
  expect_equal(length(labs2), 12)
  expect_false(any(is.na(labs2)))
  # ggplot build must not warn about missing label column
  p = suppressWarnings(plot_mv_dendrogram(r2))
  expect_s3_class(p, "ggplot")
  built = ggplot2::ggplot_build(p)
  expect_true(all(vapply(built$data, nrow, numeric(1)) > 0))
})

test_that("Davies-Bouldin uses centroid scatter (standard definition)", {
  set.seed(2)
  cl = rep(1:2, each = 15)
  X = rbind(
    cbind(rnorm(15, 0, .5), rnorm(15, 0, .5)),
    cbind(rnorm(15, 4, .5), rnorm(15, 4, .5))
  )
  dfx = as.data.frame(X)
  names(dfx) = c("a", "b")
  met = mv_cluster_metrics(dfx, a:b, cl, scale = FALSE)
  db_pkg = met$summary$value[met$summary$metric == "davies_bouldin"]
  S = vapply(1:2, function(g) {
    xg = X[cl == g, , drop = FALSE]
    cen = colMeans(xg)
    mean(sqrt(rowSums(sweep(xg, 2, cen)^2)))
  }, numeric(1))
  M = sqrt(sum((colMeans(X[cl == 1, ]) - colMeans(X[cl == 2, ]))^2))
  expect_equal(db_pkg, (S[1] + S[2]) / M, tolerance = 1e-10)
})

test_that("silhouette neighbor labels live in the original label space", {
  d5 = data.frame(
    a = c(rnorm(10), rnorm(10, 5)),
    b = c(rnorm(10), rnorm(10, 5))
  )
  cl5 = rep(c("g1", "g2"), each = 10)
  m5 = mv_cluster_metrics(d5, a:b, cl5, scale = FALSE)
  expect_true(all(m5$silhouette$neighbor %in% c("g1", "g2")))
  # every point's neighbor is the other cluster here (well-separated blobs)
  expect_true(all(m5$silhouette$neighbor != as.character(m5$silhouette$cluster)))
})

test_that("mv_corresp chi2_stat is the full-table chi-square regardless of nf", {
  set.seed(1)
  dfc = data.frame(
    a = rep(c("A", "B", "C", "D", "E"), each = 20),
    b = sample(c("u", "v", "w", "x"), 100, TRUE)
  )
  tab = table(dfc$a, dfc$b)
  chi_ref = suppressWarnings(stats::chisq.test(tab))$statistic
  r2 = mv_corresp(dfc, a, b, nf = 2)
  r3 = mv_corresp(dfc, a, b, nf = 3)
  expect_equal(r2$meta$chi2_stat, unname(chi_ref), tolerance = 1e-8)
  expect_equal(r3$meta$chi2_stat, unname(chi_ref), tolerance = 1e-8)
  expect_equal(r2$meta$chi2_stat, r3$meta$chi2_stat, tolerance = 1e-10)
})

test_that("mv_efa validates rotation", {
  expect_error(
    mv_efa(mtcars, cyl:carb, factors = 2, rotation = "varmiax"),
    "rotation"
  )
  r = mv_efa(mtcars, cyl:carb, factors = 2, rotation = "none")
  expect_s3_class(r, "mv_efa")
})

test_that("mv_kmeans validates centers and keeps center names", {
  expect_error(mv_kmeans(mtcars, cyl:carb, centers = 0), "positive integer")
  expect_error(mv_kmeans(mtcars, cyl:carb, centers = "x"), "positive integer")
  # matrix centers with wrong column count
  expect_error(
    mv_kmeans(mtcars, cyl:carb, centers = matrix(rnorm(6), 3, 2)),
    "one column per selected variable"
  )
  # named matrix centers keep row names as cluster labels (centers on the
  # raw scale to match scale = FALSE)
  raw = as.matrix(mtcars[, c("cyl", "disp", "hp")])
  cen = raw[c(1, 15, 32), , drop = FALSE]
  rownames(cen) = c("small", "mid", "big")
  r = mv_kmeans(mtcars, cyl:hp, centers = cen, scale = FALSE)
  expect_true(all(c("small", "mid", "big") %in%
    as.character(r$centers$cluster)))
  expect_true(all(as.character(r$cluster$cluster) %in%
    c("small", "mid", "big")))
  # unnamed matrix centers fall back to center_1..center_k
  cen2 = unname(raw[c(1, 15, 32), , drop = FALSE])
  stopifnot(is.null(rownames(cen2)))
  r2 = mv_kmeans(mtcars, cyl:hp, centers = cen2, scale = FALSE)
  expect_equal(as.character(r2$centers$cluster), paste0("center_", 1:3))
  # a matrix of centers must not be treated as nstart sets of centers
  expect_equal(r$meta$k, 3)
  expect_equal(sum(r$within_ss$size), 32)
})

test_that("plot_mv_elbow kmeans branch is reproducible and restores RNG", {
  set.seed(99)
  sd_before = globalenv()$.Random.seed
  g1 = plot_mv_elbow(mtcars, cyl:carb, max_k = 4)
  sd_after = globalenv()$.Random.seed
  expect_identical(sd_before, sd_after)
  g2 = plot_mv_elbow(mtcars, cyl:carb, max_k = 4)
  y1 = ggplot2::ggplot_build(g1)$data[[1]]$y
  y2 = ggplot2::ggplot_build(g2)$data[[1]]$y
  expect_equal(y1, y2)
})

test_that("mv_cluster_metrics tolerates NA labels on dropped rows only", {
  df = mtcars[, c("cyl", "disp", "hp")]
  df$cyl[1] = NA
  lab = rep(c("a", "b"), each = 16)
  lab[1] = NA # row 1 is dropped anyway (NA in cols)
  r = mv_cluster_metrics(df, cyl:hp, lab)
  expect_equal(r$meta$n_used, 31)
  # but NA among analysed rows is still an error
  lab2 = rep(c("a", "b"), each = 16)
  lab2[5] = NA
  expect_error(mv_cluster_metrics(df, cyl:hp, lab2), "NA")
})

test_that("plot_mv_dendrogram draws the k cut-height line", {
  r = mv_hclust(mtcars, cyl:carb, k = 3)
  p = plot_mv_dendrogram(r)
  built = ggplot2::ggplot_build(p)
  hline_layer = built$data[[2]]
  n = length(r$hclust$order)
  expect_equal(hline_layer$yintercept[1], r$hclust$height[n - 3])
})
