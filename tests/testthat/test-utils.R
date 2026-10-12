# Tests for the shared internal checks in R/utils.R (Phase A).

test_that(".check_matrix_like accepts data frame and matrix by default", {
  X = data.frame(a = 1:3, b = 4:6)
  expect_identical(mathmodels:::.check_matrix_like(X, "x must be a data frame or matrix."), X)
  M = as.matrix(X)
  expect_identical(mathmodels:::.check_matrix_like(M, "x must be a data frame or matrix."), M)
})

test_that(".check_matrix_like coerce = TRUE returns a data frame", {
  M = matrix(1:6, 3, 2)
  out = mathmodels:::.check_matrix_like(M, "x must be a data frame or matrix.", coerce = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(as.matrix(out), M, ignore_attr = TRUE)
})

test_that(".check_matrix_like df_only rejects a matrix with the caller's message", {
  M = matrix(1:6, 3, 2)
  expect_error(mathmodels:::.check_matrix_like(M, "data must be a data frame.", df_only = TRUE),
               "data must be a data frame", fixed = TRUE)
})

test_that(".check_matrix_like rejects non-matrix-like with the caller's message verbatim", {
  msg = "X must be a matrix or data frame."
  err = tryCatch(mathmodels:::.check_matrix_like("a", msg), error = identity)
  expect_s3_class(err, "rlang_error")
  expect_identical(conditionMessage(err), msg)
})

test_that(".check_index defaults, length and allowed-set behaviour", {
  expect_identical(mathmodels:::.check_index(NULL, 3, "L", "A"), rep(NA_character_, 3))
  expect_identical(mathmodels:::.check_index(NULL, 3, "L", "A", default = "+"), c("+", "+", "+"))
  expect_error(mathmodels:::.check_index(c("+", "-"), 3, "index must have length equal to ncol(X).", "A"),
               "must have length equal to ncol", fixed = TRUE)
  expect_error(mathmodels:::.check_index(c("+", "x"), 2, "L", 'index must contain only "+", "-", or NA'),
               'index must contain only "+", "-", or NA', fixed = TRUE)
  expect_error(mathmodels:::.check_index(c("+", NA), 2, "L", "A", allow_na = FALSE),
               "A", fixed = TRUE)
})

test_that(".check_index skips the allowed check when msg_allowed is NULL", {
  expect_identical(mathmodels:::.check_index(c("+", "-"), 2, "L", NULL), c("+", "-"))
})

test_that(".rescale_indicators is identical for matrix and data frame", {
  set.seed(1)
  M   = matrix(runif(20), 5, 4)
  idx = c("+", "-", NA, "+")
  out_m = mathmodels:::.rescale_indicators(M, idx)
  out_d = mathmodels:::.rescale_indicators(as.data.frame(M), idx)
  expect_s3_class(out_m, "data.frame")
  expect_equal(as.matrix(out_m), as.matrix(out_d), tolerance = 1e-15)
  # NA column untouched; +/- columns rescaled into [0, 1]
  expect_equal(out_d[[3]], M[, 3], tolerance = 1e-15)
  expect_true(all(out_d[[1]] >= 0 & out_d[[1]] <= 1))
  expect_true(all(out_d[[2]] >= 0 & out_d[[2]] <= 1))
})

test_that(".rescale_indicators handles all-NA index and a single column", {
  M = matrix(runif(9), 3, 3)
  out = mathmodels:::.rescale_indicators(M, c(NA, NA, NA))
  expect_equal(as.matrix(out), M, ignore_attr = TRUE, tolerance = 1e-15)
  M1 = matrix(c(3, 1, 2), 3, 1)
  out1 = mathmodels:::.rescale_indicators(M1, "+")
  expect_equal(out1[[1]], c(1, 0, 0.5), tolerance = 1e-15)
})
