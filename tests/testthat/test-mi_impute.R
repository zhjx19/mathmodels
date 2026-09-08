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
