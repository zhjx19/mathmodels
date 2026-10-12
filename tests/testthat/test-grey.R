test_that("grey_corr works", {
  ref = c(0.9, 0.8, 0.7)
  cmp = data.frame(
    x1 = c(0.9, 0.7, 0.8),
    x2 = c(0.8, 0.9, 0.7),
    x3 = c(0.7, 0.8, 0.9)
  )
  r = grey_corr(ref, cmp, rho = 0.5)
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("grey_corr with weights", {
  ref = c(0.9, 0.8, 0.7)
  cmp = data.frame(
    x1 = c(0.9, 0.7, 0.8),
    x2 = c(0.8, 0.9, 0.7),
    x3 = c(0.7, 0.8, 0.9)
  )
  r = grey_corr(ref, cmp, rho = 0.5, w = c(0.3, 0.3, 0.4))
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("grey_corr_input validation", {
  ref = c(0.9, 0.8)
  cmp = data.frame(x1 = c(0.9, 0.7), x2 = c(0.8, 0.9))
  expect_error(grey_corr("a", cmp), "must be a numeric vector")
  expect_error(grey_corr(ref, cmp, w = c(0.5)), "must have length")
})

test_that("grey_corr_topsis works", {
  cmp = data.frame(
    x1 = c(8, 7, 6),
    x2 = c(150, 180, 200),
    x3 = c(60, 80, 100)
  )
  w = c(0.3, 0.4, 0.3)
  idx = c("+", "+", "+")
  r = grey_corr_topsis(cmp, w, idx, rho = 0.5)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("grey_corr_topsis handles non-square matrices (n != m)", {
  # 4 alternatives x 3 indicators: regression for the lost transpose, which
  # errored with "Length of ref must equal number of rows in cmp."
  X = data.frame(
    x1 = c(8, 7, 6, 9),
    x2 = c(150, 180, 200, 160),
    x3 = c(60, 80, 100, 70)
  )
  w   = c(0.3, 0.4, 0.3)
  idx = c("+", "+", "+")
  r = grey_corr_topsis(X, w, idx, rho = 0.5)
  expect_length(r, nrow(X))
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(r, c(0.3589743590, 0.4901108455, 0.5806451613, 0.5018455581),
               tolerance = 1e-8)
})

test_that("grey_corr_topsis returns one degree per alternative for any shape", {
  set.seed(1)
  for (d in list(c(4, 3), c(20, 4), c(2, 5), c(3, 3))) {
    X = as.data.frame(matrix(runif(d[1] * d[2]), d[1], d[2]))
    r = grey_corr_topsis(X, rep(1 / d[2], d[2]), rep("+", d[2]))
    expect_length(r, d[1])
    expect_true(all(r >= 0 & r <= 1))
  }
})

test_that("grey_corr_topsis accepts matrix input like a data frame", {
  set.seed(1)
  M   = matrix(runif(60), 12, 5)
  w   = c(0.1, 0.25, 0.2, 0.15, 0.3)
  idx = c("+", "-", NA, "+", "-")
  expect_equal(grey_corr_topsis(M, w, idx),
               grey_corr_topsis(as.data.frame(M), w, idx),
               tolerance = 1e-12)
})

test_that("grey_corr_topsis skips rescaling when every index is NA", {
  X = data.frame(x1 = c(3, 1, 2), x2 = c(9, 5, 7), x3 = c(2, 4, 6))
  r = grey_corr_topsis(X, c(1, 1, 1), NULL)
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("grey_corr_topsis input validation", {
  cmp = data.frame(x1=c(1,2,4), x2=c(2,3,5))
  expect_error(grey_corr_topsis(cmp, c(0.5, 0.3, 0.2), NULL), "must have length")
})

test_that("grey_corr_topsis enforces the same size floors as topsis", {
  expect_error(grey_corr_topsis(data.frame(x = 1:3), 1), "at least 2 columns")
  expect_error(grey_corr_topsis(data.frame(a = 1, b = 2), c(0.5, 0.5)), "at least 2 rows")
})

test_that("grey_corr_topsis reports readable errors for NA and non-finite w", {
  X = data.frame(x1 = c(8, 7, 6, 9), x2 = c(150, 180, 200, 160), x3 = c(60, 80, 100, 70))
  expect_error(grey_corr_topsis(X, c(0.4, 0.3, NA)), "no NA")
  expect_error(grey_corr_topsis(X, c(0.4, 0.3, Inf)), "finite")
})
