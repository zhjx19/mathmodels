# Tests for the FCS engine and impute_multiple() (mi_impute.R)

test_that(".fcs_methods assigns by column type", {
  d = data.frame(
    x = c(1, NA, 3, 4, 5),
    z = c(1.5, 2, NA, 4, 5),
    b = factor(c("y", NA, "n", "y", "n")),     # 2 levels
    p = factor(c("a", "b", NA, "c", "a")),      # 3 levels
    k = c("u", NA, "v", "u", "v")               # character, 2 levels
  )
  mth = mathmodels:::.fcs_methods(d, names(d))
  expect_equal(unname(mth),
               c("pmm", "pmm", "logreg", "polyreg", "logreg"))
  expect_equal(names(mth), names(d))
})

test_that(".fcs_init fills numeric with mean, categorical with mode", {
  d = data.frame(x = c(1, NA, 3), g = factor(c("a", NA, "a")))
  init = mathmodels:::.fcs_init(d, c("x", "g"))
  expect_false(anyNA(init$x))
  expect_equal(init$x[2], 2)
  expect_equal(as.character(init$g[2]), "a")
  # original untouched
  expect_true(is.na(d$x[2]))
})

test_that(".fcs_design dummy-encodes factors, keeps intercept, drops NA features", {
  d = data.frame(
    y = c(1, NA, 3, 4),
    x = c(2, 3, NA, 5),                          # has NA -> not usable yet
    g = factor(c("a", "b", "a", "b"))
  )
  des = mathmodels:::.fcs_design(d, d, "y")
  expect_equal(des$na_idx, 2)
  expect_equal(des$obs_idx, c(1, 3, 4))
  # X: intercept + g dummy (x excluded because df_cur still has NA there)
  expect_equal(ncol(des$X), 2)
  expect_type(des$y_obs, "double")
})

test_that(".fcs_bayes_lm draws beta around OLS fit with positive sigma2", {
  set.seed(1)
  X = cbind(1, 1:30)
  y = 2 + 0.5 * (1:30) + rnorm(30)
  draw = mathmodels:::.fcs_bayes_lm(X, y)
  expect_length(draw$beta, 2)
  expect_gt(draw$sigma2, 0)
  # repeated draws average near the OLS coefficients (Bayes draws unbiased)
  draws = replicate(2000, mathmodels:::.fcs_bayes_lm(X, y)$beta)
  expect_equal(rowMeans(draws), unname(stats::lm.fit(X, y)$coefficients),
               tolerance = 0.05)
})

test_that(".fcs_impute_col pmm keeps imputed values inside observed range", {
  set.seed(1)
  d = data.frame(y = c(1, 2, 3, 10, NA), x = c(1, 2, 3, 4, 4.5))
  # donors are observed y values -> the mean of 3 closest never leaves [1, 10]
  out = mathmodels:::.fcs_impute_col(
    d, d, "y", "pmm", new.env(parent = emptyenv()))
  expect_false(anyNA(out))
  expect_gte(min(out), 1)
  expect_lte(max(out), 10)
})

test_that(".fcs_impute_col norm adds noise and stays finite", {
  set.seed(1)
  d = data.frame(y = c(1, 2, 3, 10, NA, 4), x = c(1, 2, 3, 4, 4.5, 5))
  out = mathmodels:::.fcs_impute_col(
    d, d, "y", "norm", new.env(parent = emptyenv()))
  expect_false(anyNA(out))
  expect_true(all(is.finite(out)))
})

test_that(".fcs_impute_col logreg returns a binary factor with original levels", {
  set.seed(1)
  d = data.frame(
    g = factor(c("y", "n", "y", "n", NA, "y")),
    x = c(1, 2, 1.5, 2.2, 1.2, 0.8)
  )
  out = mathmodels:::.fcs_impute_col(
    d, d, "g", "logreg", new.env(parent = emptyenv()))
  expect_s3_class(out, "factor")
  expect_setequal(levels(out), c("n", "y"))
  expect_false(anyNA(out))
})

test_that(".fcs_impute_col polyreg keeps all original levels available", {
  set.seed(1)
  d = data.frame(
    g = factor(c("a", "b", "c", "a", NA, "b", "c", "a")),
    x = c(1, 2, 3, 1.5, 2.5, 2.2, 3.1, 0.9)
  )
  out = mathmodels:::.fcs_impute_col(
    d, d, "g", "polyreg", new.env(parent = emptyenv()))
  expect_s3_class(out, "factor")
  expect_setequal(levels(out), c("a", "b", "c"))
  expect_false(anyNA(out))
})

# ---------------------------------------------------------------------------
# impute_multiple
# ---------------------------------------------------------------------------

make_mi_data = function(n = 60) {
  set.seed(7)
  x = rnorm(n)
  z = 0.5 * x + rnorm(n)
  g = factor(ifelse(x + rnorm(n) > 0, "a", "b"))
  d = data.frame(x = x, z = z, g = g)
  d$x[sample(n, 8)] = NA
  d$z[sample(n, 6)] = NA
  d$g[sample(n, 5)] = NA
  d
}

test_that("impute_multiple: shape, types, and factor levels preserved", {
  d = make_mi_data()
  imps = impute_multiple(d, m = 3, maxit = 3, seed = 42)
  expect_s3_class(imps, "mi_impute")
  expect_length(imps, 3)
  for (i in 1:3) {
    expect_equal(dim(imps[[i]]), dim(d))
    expect_setequal(names(imps[[i]]), names(d))
    expect_s3_class(imps[[i]]$g, "factor")
    expect_setequal(levels(imps[[i]]$g), c("a", "b"))
    expect_false(anyNA(imps[[i]]))
  }
})

test_that("impute_multiple: imputed cells only touch originally-missing rows", {
  d = make_mi_data()
  imps = impute_multiple(d, m = 2, maxit = 3, seed = 1)
  for (i in 1:2) {
    expect_equal(imps[[i]]$x[!is.na(d$x)], d$x[!is.na(d$x)])
    expect_equal(imps[[i]]$z[!is.na(d$z)], d$z[!is.na(d$z)])
    expect_equal(imps[[i]]$g[!is.na(d$g)], d$g[!is.na(d$g)])
  }
})

test_that("impute_multiple: logical column with NA comes back as factor", {
  d = data.frame(
    b = c(TRUE, FALSE, NA, TRUE, FALSE, NA, TRUE, FALSE, TRUE, FALSE,
          NA, TRUE, FALSE, TRUE, FALSE, NA, TRUE, FALSE, TRUE, FALSE),
    x = c(1.2, 0.3, 0.8, 2.1, 0.5, 1.1, 1.8, 0.2, 0.9, 1.5,
          2.0, 0.4, 1.3, 0.7, 1.9, 0.6, 1.4, 0.1, 1.6, 1.0)
  )
  imps = impute_multiple(d, m = 2, maxit = 2, seed = 11)
  expect_s3_class(imps, "mi_impute")
  expect_setequal(levels(imps[[1]]$b), c("FALSE", "TRUE"))
  expect_false(anyNA(imps[[1]]$b))
  expect_false(anyNA(imps[[2]]$b))
})

test_that("impute_multiple: seed reproducibility and between-set variety", {
  d = make_mi_data()
  a = impute_multiple(d, m = 2, maxit = 3, seed = 99)
  b = impute_multiple(d, m = 2, maxit = 3, seed = 99)
  expect_equal(a[[1]]$x, b[[1]]$x)
  expect_equal(a[[2]]$g, b[[2]]$g)
  # two imputations differ on at least some numeric cells
  expect_false(isTRUE(all.equal(a[[1]]$x[is.na(d$x)],
                                a[[2]]$x[is.na(d$x)])))
})

test_that("impute_multiple: attributes carry methods, trace, n_imputed", {
  d = make_mi_data()
  imps = impute_multiple(d, m = 2, maxit = 2, seed = 5)
  mth = attr(imps, "methods")
  expect_true(all(mth == c(x = "pmm", z = "pmm", g = "logreg")))
  expect_equal(attr(imps, "iterations"), 2L)
  tr = attr(imps, "trace")
  expect_s3_class(tr, "tbl_df")
  expect_equal(names(tr),
               c("imputation", "iteration", "variable", "mean", "sd"))
  expect_equal(nrow(tr), 2 * 2 * 3)  # m * maxit * n_imputed_cols
  expect_equal(attr(imps, "n_imputed"), c(x = 8L, z = 6L, g = 5L))
})

test_that("impute_multiple: methods override and .cols selection", {
  d = make_mi_data()
  imps = impute_multiple(d, methods = c(x = "norm"), m = 2, maxit = 2,
                         seed = 3)
  expect_equal(attr(imps, "methods")[["x"]], "norm")
  # .cols limits imputation to the selected column; z keeps its NAs
  imps2 = impute_multiple(d, .cols = x, m = 2, maxit = 2, seed = 3)
  expect_true(anyNA(imps2[[1]]$z))
  expect_false(anyNA(imps2[[1]]$x))
  expect_equal(names(attr(imps2, "methods")), "x")
})

test_that("impute_multiple: input validation", {
  d = make_mi_data()
  expect_error(impute_multiple(d, m = 1), "m must be")
  expect_error(impute_multiple(d, maxit = 0), "maxit must be")
  expect_error(impute_multiple(d, methods = c(x = "rf")), "unknown")
  expect_error(impute_multiple(d, .cols = y), "no existing")
  expect_error(impute_multiple(as.matrix(d)), "data.frame")
})

test_that("impute_multiple: no-missing data returns copies unchanged", {
  d = data.frame(x = 1:5, g = factor(c("a", "b", "a", "a", "b")))
  imps = impute_multiple(d, m = 2, maxit = 2, seed = 1)
  expect_equal(imps[[1]], d)
  expect_equal(attr(imps, "n_imputed"), setNames(integer(), character()))
})

test_that("impute_multiple: print method runs", {
  d = make_mi_data()
  imps = impute_multiple(d, m = 2, maxit = 2, seed = 1)
  out = capture.output(print(imps))
  expect_true(any(grepl("2 imputation", out)))
})

test_that("impute_multiple norm agrees with mice norm in pooled inference", {
  skip_if_not_installed("mice")
  set.seed(11)
  n = 120
  x = rnorm(n)
  z = 0.7 * x + rnorm(n)
  d = data.frame(x = x, z = z)
  d$y = 0.5 * x + 0.3 * z + rnorm(n)
  d$x[sample(n, 20)] = NA
  d$z[sample(n, 15)] = NA

  imps = impute_multiple(d, .cols = c(x, z), m = 10, maxit = 5, seed = 3,
                         methods = c(x = "norm", z = "norm"))
  fits = lapply(imps, \(dat) lm(y ~ x + z, dat))
  ours = pool_fit(fits)

  mids = mice::mice(d, m = 10, maxit = 5, seed = 3, printFlag = FALSE,
                    method = c(x = "norm", z = "norm", y = ""))
  fits_mice = lapply(seq_len(10), \(i) lm(y ~ x + z, mice::complete(mids, i)))
  ref = mice::pool(fits_mice)

  cmp = merge(ours[, c("term", "estimate", "std_error")],
              summary(ref)[, c("term", "estimate", "std.error")],
              by = "term")
  for (trm in cmp$term) {
    row = cmp[cmp$term == trm, ]
    expect_equal(row$estimate.x, row$estimate.y, tolerance = 0.15)
    expect_equal(row$std_error, row$std.error, tolerance = 0.15)
  }
})
