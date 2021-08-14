test_that("AHP weights and type", {
  A = matrix(c(1,   1/2,
               2,   1), byrow = TRUE, nrow = 2)
  rlt = AHP(A)
  expect_equal(rlt$W, c(0.3333, 0.6667), tolerance = 0.001)
  expect_type(rlt, "list")
})
