test_that("weight_entropy works", {
  X = data.frame(
    x1 = c(3, 5, 2, 7),
    x2 = c(10, 20, 15, 25)
  )
  index = c("+", "-")
  r = weight_entropy(X, index)
  expect_type(r, "list")
  expect_named(r, c("w", "s"))
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 4)
})

test_that("weight_entropy defaults", {
  X = data.frame(x1=c(3,5,2,7), x2=c(0.2,0.3,0.4,0.5))
  r = weight_entropy(X)
  expect_equal(sum(r$w), 1, tolerance=1e-8)
})

test_that("weight_entropy input validation", {
  X = data.frame(x1 = c(3, 5, 2, 7))
  expect_error(weight_entropy(X), "at least 2 columns")
  expect_error(weight_entropy("not a df"), "must be a data frame or matrix")
})

test_that("weight_entropy matrix input equals data.frame input", {
  Xdf = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  index = c("+", "-")
  r1 = weight_entropy(Xdf, index)
  r2 = weight_entropy(as.matrix(Xdf), index)
  expect_equal(r1$w, r2$w)
  expect_equal(r1$s, r2$s)
})

test_that("weight_entropy negative and NA data errors", {
  expect_error(weight_entropy(data.frame(x1=c(3,-5,2,7), x2=c(1,2,3,4))), "non-negative")
  expect_error(weight_entropy(data.frame(x1=c(3,5,2,NA), x2=c(1,2,3,4))), "NA")
})

test_that("weight_entropy index validation", {
  X = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  expect_error(weight_entropy(X, c("+", "pos")), 'index must contain only')
})

test_that("weight_entropy warns when NA and rescaled columns mix scales", {
  X = data.frame(a = 1:10, b = 1:10)   # a is "+" (rescaled), b is NA (raw)
  expect_warning(weight_entropy(X, index = c("+", NA)), "already normalized")
})

test_that("weight_entropy stays quiet when NA columns are normalised", {
  X = data.frame(a = 1:10, b = seq(0.1, 1, length.out = 10))
  expect_no_warning(weight_entropy(X, index = c("+", NA)))
  # + and - columns are rescaled internally, so they never trip the guard
  Y = data.frame(x1 = c(3, 5, 2, 7), x2 = c(10, 20, 15, 25))
  expect_no_warning(weight_entropy(Y, c("+", "-")))
  # index = NULL rescales nothing, so raw data is not reported
  expect_no_warning(weight_entropy(data.frame(x1 = c(3, 5, 2, 7),
                                              x2 = c(10, 20, 15, 25))))
})

