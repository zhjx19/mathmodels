test_that("linear_sum works", {
  data = data.frame(a=c(0.8,0.7,0.9), b=c(0.6,0.8,0.5))
  w = c(0.4, 0.6)
  r = linear_sum(data, w)
  expect_length(r, 3)
  expect_equal(r[1], 0.8*0.4 + 0.6*0.6)
})

test_that("linear_sum input validation", {
  data = data.frame(a=c(0.8,0.7), b=c(0.6,0.8))
  expect_error(linear_sum(data, c(0.4, 0.3, 0.3)), "must equal ncol")
  expect_error(linear_sum("a", c(1,2)), "must be a data frame or matrix")
})
