test_that("entropy_weight works", {
  X = data.frame(
    x1 = c(3, 5, 2, 7),
    x2 = c(10, 20, 15, 25)
  )
  index = c("+", "-")
  r = entropy_weight(X, index)
  expect_type(r, "list")
  expect_named(r, c("w", "s"))
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 4)
})

test_that("entropy_weight defaults", {
  X = data.frame(x1=c(3,5,2,7), x2=c(0.2,0.3,0.4,0.5))
  r = entropy_weight(X)
  expect_equal(sum(r$w), 1, tolerance=1e-8)
})

test_that("entropy_weight input validation", {
  X = data.frame(x1 = c(3, 5, 2, 7))
  expect_error(entropy_weight(X), "at least 2 columns")
  expect_error(entropy_weight("not a df"), "must be a data frame or matrix")
})

test_that("entropy_weight matrix input equals data.frame input", {
  Xdf = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  index = c("+", "-")
  r1 = entropy_weight(Xdf, index)
  r2 = entropy_weight(as.matrix(Xdf), index)
  expect_equal(r1$w, r2$w)
  expect_equal(r1$s, r2$s)
})

test_that("entropy_weight negative and NA data errors", {
  expect_error(entropy_weight(data.frame(x1=c(3,-5,2,7), x2=c(1,2,3,4))), "non-negative")
  expect_error(entropy_weight(data.frame(x1=c(3,5,2,NA), x2=c(1,2,3,4))), "NA")
})

test_that("entropy_weight index validation", {
  X = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  expect_error(entropy_weight(X, c("+", "pos")), 'index must contain only')
})

