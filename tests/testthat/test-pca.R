test_that("weight_pca works", {
  ind = c("+","+","-","-")
  r = weight_pca(iris[1:10, 1:4], ind, nfs = 2)
  expect_type(r, "list")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 10)
})

test_that("weight_pca input validation", {
  X = data.frame(x1=c(3,5,2,7))
  expect_error(weight_pca(X), "at least 2 columns")
})

test_that("weight_pca treats matrix and data frame input identically", {
  set.seed(1)
  M   = matrix(runif(120), 30, 4)
  ind = c("+", "+", "-", "-")
  r_df = weight_pca(as.data.frame(M), ind, nfs = 2)
  r_m  = weight_pca(M, ind, nfs = 2)
  expect_equal(r_m$w, r_df$w, tolerance = 1e-12)
  expect_equal(r_m$s, r_df$s, tolerance = 1e-12)
})

test_that("weight_fa matches factanal loadings-based weights", {
  set.seed(1)
  X = data.frame(
    x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50),
    x4 = rnorm(50), x5 = rnorm(50)
  )
  r = weight_fa(X, nfs = 2)
  expect_type(r, "list")
  expect_equal(sum(r$w), 1, tolerance = 1e-10)
  # hand-computed reference: abs(loadings) %*% variance shares, normalized
  fit = stats::factanal(X, factors = 2, rotation = "varimax", scores = "none")
  A = unclass(fit$loadings)
  share = colSums(A^2) / 5
  w_ref = as.vector(abs(A) %*% share)
  w_ref = w_ref / sum(w_ref)
  expect_equal(r$w, w_ref, tolerance = 1e-12)
  expect_equal(r$share, as.vector(share), tolerance = 1e-12)
  expect_equal(r$uniqueness, as.numeric(fit$uniquenesses), tolerance = 1e-12)
  # scores = X %*% w
  expect_equal(r$s, as.vector(as.matrix(X) %*% w_ref), tolerance = 1e-10)
  # loading matrix shape
  expect_equal(dim(r$loading), c(5, 2))
})

test_that("weight_fa direction index flips signs only", {
  set.seed(2)
  X = data.frame(a = rnorm(40), b = rnorm(40), c = rnorm(40))
  r0 = weight_fa(data.frame(
    a = scale(X$a), b = -scale(X$b), c = scale(X$c)), nfs = 1)
  r1 = weight_fa(X, index = c("+", "-", NA), nfs = 1)
  # with abs loadings, flipping a column sign must not change the weights
  expect_equal(as.vector(r1$w), as.vector(r0$w), tolerance = 1e-10)
  expect_equal(sum(r1$w), 1, tolerance = 1e-10)
})

test_that("weight_fa input validation", {
  set.seed(3)
  X = data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  expect_error(weight_fa(X, nfs = 3), "smaller")
  expect_error(weight_fa(X[1:3, ], nfs = 1), "more rows")
  expect_error(weight_fa(X, nfs = 0), "positive integer")
  expect_error(weight_fa(X, rotation = "varmiax"), "rotation")
  expect_error(weight_fa(X, method = "cube"), 'method must be "abs"')
  Xna = X; Xna$a[1] = NA
  expect_error(weight_fa(Xna), "NA")
  Xconst = X; Xconst$c = 5
  expect_error(weight_fa(Xconst), "constant columns")
})

test_that("weight_fa index validation keeps its exact messages", {
  # Pinning test: these two assertions pass BEFORE the migration to
  # .check_index() and must still pass AFTER it -- the task is de-duplication,
  # not a behaviour change.
  set.seed(1)
  X = data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  expect_error(weight_fa(X, index = c("+", "+")),
               "index must have length equal to ncol(X), or be NULL.", fixed = TRUE)
  expect_error(weight_fa(X, index = c("+", "x", "+")),
               'index must contain only "+", "-", or NA', fixed = TRUE)
})
