test_that("cv_weight works", {
  X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
  r = cv_weight(X)
  expect_equal(sum(r), 1, tolerance=1e-8)
  expect_length(r, 2)
})

test_that("cv_weight input validation", {
  expect_error(cv_weight(data.frame(x=c(10,20))), "at least 2 columns")
  expect_error(cv_weight(data.frame(x=c(0,10), y=c(5,20))), "must be positive")
})
