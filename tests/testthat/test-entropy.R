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
