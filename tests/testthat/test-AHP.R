test_that("AHP weights and type", {
  A = matrix(c(1,   1/2,
               2,   1), byrow = TRUE, nrow = 2)
  rlt = AHP(A)
  expect_equal(rlt$w, c(0.3333, 0.6667), tolerance = 0.001)
  expect_type(rlt, "list")
})

test_that("AHP consistency", {
  A = matrix(c(1,   1/2, 4, 3,   3,
               2,   1,   7, 5,   5,
               1/4, 1/7, 1, 1/2, 1/3,
               1/3, 1/5, 2, 1,   1,
               1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
  rlt = AHP(A)
  expect_equal(rlt$CR, rlt$CI / c(0,0,0.58,0.90,1.12,1.24,1.32,1.41,1.45,1.49,
                                   1.51,1.48,1.56,1.57,1.59)[5])
  expect_true(rlt$Lmax > 0)
  expect_equal(sum(rlt$w), 1, tolerance = 1e-8)
})

test_that("AHP input validation", {
  expect_error(AHP("not a matrix"), "must be a matrix")
  expect_error(AHP(matrix(1:6, nrow=2)), "must be a square matrix")
  expect_error(AHP(matrix(1, nrow=16, ncol=16)), "at most 15")
})
