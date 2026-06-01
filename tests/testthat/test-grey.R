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

test_that("grey_corr_topsis input validation", {
  cmp = data.frame(x1=c(1,2,4), x2=c(2,3,5))
  expect_error(grey_corr_topsis(cmp, c(0.5, 0.3, 0.2), NULL), "must have length")
})
