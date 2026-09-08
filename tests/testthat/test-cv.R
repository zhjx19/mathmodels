test_that("weight_cv works", {
  X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
  r = weight_cv(X)
  expect_equal(sum(r), 1, tolerance=1e-8)
  expect_length(r, 2)
})

test_that("weight_cv input validation", {
  expect_error(weight_cv(data.frame(x=c(10,20))), "at least 2 columns")
  expect_error(weight_cv(data.frame(x=c(0,10), y=c(5,20))), "must be positive")
  expect_error(weight_cv(data.frame(x=c(NA,10), y=c(5,20))), "NA")
})

