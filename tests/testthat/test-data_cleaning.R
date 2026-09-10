# Tests for the data cleaning toolkit (data_cleaning.R)

set.seed(42)
d_ok = data.frame(g = rep(c("a", "b"), each = 5),
                  x = c(1, 2, NA, 4, 5, 10, NA, 12, 14, 16),
                  y = c(2, 4, 6, 8, 10, 1, 3, NA, 7, 9))

# ---------------------------------------------------------------------------
# na_summary
# ---------------------------------------------------------------------------

test_that("na_summary counts per column", {
  r = na_summary(d_ok, .cols = x:y)
  expect_equal(r$variable, c("x", "y"))
  expect_equal(r$n_missing, c(2, 1))
  expect_equal(r$pct_missing, c(0.2, 0.1), tolerance = 1e-12)
  expect_equal(nrow(na_summary(d_ok)), 3)  # default: all columns
})

# ---------------------------------------------------------------------------
# na_mcar_test
# ---------------------------------------------------------------------------

test_that("na_mcar_test: random missing gives large p, correlated gives small", {
  set.seed(1)
  d_rand = data.frame(x = rnorm(200), y = rnorm(200), z = rnorm(200))
  d_rand$x[sample(200, 20)] = NA
  r = na_mcar_test(d_rand, .cols = x:z)
  expect_s3_class(r, "tbl_df")
  expect_equal(r$df, 3 + 2 - 3)  # sum(pattern-observed cols) - p
  expect_gt(r$p_value, 0.01)

  # strongly non-random: y missing exactly on the 10 planted extreme
  # values of x -> the pattern's mean vector is far from the overall mean
  set.seed(1)
  d_dep = data.frame(x = c(rnorm(190), rep(30, 10)),
                     y = rnorm(200), z = rnorm(200))
  d_dep$y[d_dep$x > 15] = NA
  r2 = na_mcar_test(d_dep, .cols = x:z)
  expect_lt(r2$p_value, 0.001)
  # patterns attribute
  expect_true(nrow(attr(r, "patterns")) >= 1)
})

# plain correlated numeric matrix
rmvlike = function(n = 80) {
  x = rnorm(n); y = 2 * x + rnorm(n); z = -x + y + rnorm(n)
  cbind(x, y, z)
}

test_that("na_mcar_test EM estimate matches complete-case moments", {
  set.seed(3)
  X = rmvlike(80)
  em = mathmodels:::.em_norm_mle(X)
  expect_equal(em$mu, colMeans(X), tolerance = 1e-8)
  expect_equal(em$sigma, stats::cov(X) * 79 / 80, tolerance = 1e-6)
})

test_that("na_mcar_test validation", {
  expect_error(na_mcar_test(data.frame(x = 1:5), .cols = x), "no missing")
  expect_warning(
    na_mcar_test(data.frame(x = c(NA_real_, NA_real_), y = 1:2), .cols = x:y),
    "fully missing")
})

# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

test_that("na plots return ggplot with data", {
  p1 = plot_na_heatmap(d_ok, .cols = x:y)
  p2 = plot_na_bar(d_ok, .cols = x:y)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_true(any(vapply(ggplot2::ggplot_build(p1)$data, nrow, numeric(1)) > 0))
  expect_true(any(vapply(ggplot2::ggplot_build(p2)$data, nrow, numeric(1)) > 0))
})

# ---------------------------------------------------------------------------
# impute
# ---------------------------------------------------------------------------

test_that("impute mean/median fill correctly and record changes", {
  r = impute(d_ok, .cols = x, method = "mean")
  expect_equal(r$x[3], mean(d_ok$x, na.rm = TRUE), tolerance = 1e-12)
  expect_equal(sum(is.na(r$x)), 0)
  ch = attr(r, "changes")
  expect_equal(nrow(ch), 2)
  expect_equal(sort(ch$row), c(3, 7))
  expect_true(all(is.na(ch$old)))
  r2 = impute(d_ok, .cols = x:y, method = "median")
  expect_equal(r2$y[8], stats::median(d_ok$y, na.rm = TRUE), tolerance = 1e-12)
  expect_equal(nrow(attr(r2, "changes")), 3)
})

test_that("impute grouped mean uses .by slices", {
  r = impute(d_ok, .cols = x, method = "mean", .by = g)
  xa = d_ok$x[d_ok$g == "a"]
  expect_equal(r$x[3], mean(xa, na.rm = TRUE), tolerance = 1e-12)
  xb = d_ok$x[d_ok$g == "b"]
  expect_equal(r$x[7], mean(xb, na.rm = TRUE), tolerance = 1e-12)
})

test_that("impute mode works for numeric and character columns", {
  d = data.frame(num = c(1, 2, 2, NA), txt = c("a", "b", "a", NA))
  r = impute(d, .cols = num:txt, method = "mode")
  expect_equal(r$num[4], 2)
  expect_equal(as.character(r$txt[4]), "a")
  # numeric column stays numeric
  expect_true(is.numeric(r$num))
})

test_that("impute linear interpolation recovers straight lines exactly", {
  d = data.frame(t = 1:10, x = (1:10) * 3 + 1)
  d$x[c(4, 7)] = NA
  r = impute(d, .cols = x, method = "linear")
  expect_equal(r$x, (1:10) * 3 + 1, tolerance = 1e-10)
  # end values cannot be interpolated
  d2 = data.frame(x = c(NA, 1, 2, 3))
  r2 = impute(d2, .cols = x, method = "linear")
  expect_true(is.na(r2$x[1]))
})

test_that("impute constant and validation", {
  r = impute(d_ok, .cols = x, method = "constant", value = 0)
  expect_equal(r$x[3], 0)
  expect_error(impute(d_ok, .cols = x, method = "constant"), "value")
  expect_error(impute(d_ok, .cols = x, method = "linear", .by = g),
    ".by grouping")
  dfc = data.frame(x = c("a", NA, "b"))
  expect_error(impute(dfc, .cols = x, method = "mean"), "numeric")
})

test_that("impute leaves untouched columns and no-missing data as is", {
  r = impute(d_ok, .cols = x, method = "mean")
  expect_equal(r$g, d_ok$g)
  expect_equal(r$y, d_ok$y)
  d_no_na = data.frame(x = 1:5)
  r2 = impute(d_no_na, .cols = x, method = "mean")
  expect_equal(r2$x, 1:5)
  expect_equal(nrow(attr(r2, "changes")), 0)
})

# ---------------------------------------------------------------------------
# impute_model
# ---------------------------------------------------------------------------

test_that("impute_model lm recovers linear structure", {
  set.seed(1)
  d = data.frame(x1 = 1:40, x2 = rnorm(40))
  d$y = 2 * d$x1 + 3 * d$x2 + rnorm(40, sd = 0.1)
  d$y[c(5, 20)] = NA
  r = impute_model(d, .cols = y, model = "lm")
  expect_equal(sum(is.na(r$y)), 0)
  # recovered values close to the true linear structure
  true_vals = 2 * d$x1[c(5, 20)] + 3 * d$x2[c(5, 20)]
  expect_equal(r$y[c(5, 20)], true_vals, tolerance = 0.5)
  det = attr(r, "details")
  expect_equal(det$n_train, 38)
  expect_equal(det$n_imputed, 2)
})

test_that("impute_model knn fills from neighbors", {
  set.seed(2)
  d = data.frame(
    f1 = c(rnorm(20, 0), rnorm(20, 10)),
    f2 = c(rnorm(20, 0), rnorm(20, 10))
  )
  d$y = c(rnorm(20, 0), rnorm(20, 10))
  d$y[40] = NA  # a point in the second cluster
  r = impute_model(d, .cols = y, model = "knn", k = 3)
  expect_equal(sum(is.na(r$y)), 0)
  expect_gt(r$y[40], 5)  # filled from the high cluster
  det = attr(r, "details")
  expect_equal(det$model, "knn")
})

test_that("impute_model tree runs and fills", {
  set.seed(3)
  d = data.frame(x1 = rnorm(30), x2 = rnorm(30))
  d$y = d$x1 + rnorm(30, sd = 0.2)
  d$y[7] = NA
  r = impute_model(d, .cols = y, model = "tree")
  expect_false(is.na(r$y[7]))
  expect_equal(sum(is.na(r$y)), 0)
})

test_that("impute_model .by slices and validation", {
  set.seed(4)
  d = data.frame(
    g = rep(c("a", "b"), each = 15),
    x1 = c(rnorm(15, 0), rnorm(15, 100)),
    y = c(rnorm(15, 0), rnorm(15, 100))
  )
  d$y[c(3, 18)] = NA
  r = impute_model(d, .cols = y, model = "knn", .by = g, k = 3)
  expect_equal(sum(is.na(r$y)), 0)
  expect_lt(abs(r$y[3]), 20)   # filled within the low cluster
  expect_gt(r$y[18], 80)       # filled within the high cluster
  expect_error(impute_model(d, .cols = y, model = "knn", k = 0), "positive")
  d_no_feat = data.frame(y = c(1, NA, 3))
  expect_error(impute_model(d_no_feat, .cols = y, model = "lm"),
    "no numeric feature")
})

# ---------------------------------------------------------------------------
# outlier
# ---------------------------------------------------------------------------

test_that("outlier detect does not touch the data and reports details", {
  d = data.frame(g = rep(c("a", "b"), each = 20),
                 x = c(rnorm(20, 0), rnorm(19, 10), 100))
  r = outlier(d, .cols = x, .by = g)  # default: detect
  expect_s3_class(r, "tbl_df")        # detail table, not the data
  expect_true(all(c("variable", "row", "value", "bound_low",
    "bound_high") %in% names(r)))
  hit = r$row[r$variable == "x"]
  expect_true(100 %in% d$x[hit])  # the planted outlier is among the hits
  # summary attribute
  summ = attr(r, "summary")
  expect_true(all(c("variable", "n_outliers", "method") %in% names(summ)))
})

test_that("outlier to='na' sets NA and to='cap' caps at fence", {
  d = data.frame(x = c(rnorm(19, 10), 100))
  r_na = outlier(d, .cols = x, to = "na")
  expect_true(is.na(r_na$x[20]))
  det = attr(r_na, "outliers")
  expect_equal(det$row, 20)
  # cap: value equals the upper fence
  r_cap = outlier(d, .cols = x, to = "cap")
  expect_equal(r_cap$x[20], det$bound_high, tolerance = 1e-12)
  expect_false(is.na(r_cap$x[20]))
})

test_that("outlier .by computes bounds within slices", {
  d = data.frame(
    g = rep(c("a", "b"), each = 15),
    x = c(rnorm(14, 0), 50, rnorm(14, 100), 50)
  )
  r = outlier(d, .cols = x, .by = g)
  # two outliers: one per slice (50 vs the opposite cluster)
  expect_equal(nrow(r), 2)
  a_hit = r$row[r$g == "a"]
  expect_equal(d$x[a_hit], 50)
  expect_gt(r$bound_low[r$g == "b"], 40)
})

test_that("outlier constant columns are skipped with warning", {
  d = data.frame(x = rep(5, 10), y = c(rnorm(9), 100))
  o = withCallingHandlers(
    outlier(d, .cols = x:y),
    warning = function(w) {
      expect_match(conditionMessage(w), "constant")
      invokeRestart("muffleWarning")
    })
  expect_equal(nrow(o), 1)  # only y's outlier
})

test_that("outlier validation", {
  expect_error(outlier(d_ok, .cols = x, k = -1), "positive")
  expect_error(outlier(data.frame(x = c("a", "b")), .cols = x), "numeric")
})
