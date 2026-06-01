test_that("topsis works", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  w = c(0.6, 0.4)
  idx = c("+","-")
  r = topsis(A, w, idx)
  expect_true(all(r >= 0 & r <= 1))
  expect_length(r, 3)
})

test_that("topsis defaults", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  r = topsis(A)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("topsis input validation", {
  expect_error(topsis(data.frame(x=c(1,2))), "at least 2 columns")
  expect_error(topsis(data.frame(x1=c(2,5,3), x2=c(8,1,6)), w=c(0.5)), "must have length equal to ncol")
})
