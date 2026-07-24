# =============================================================================
# Unit tests for the interpolation & curve fitting module interp_fit.R
# =============================================================================

# --- ._validate_xy ---

test_that("._validate_xy() passes for valid input", {
  expect_no_error(mathmodels:::`._validate_xy`(c(1,2,3), c(4,5,6)))
})

test_that("._validate_xy() rejects non-numeric x", {
  expect_error(mathmodels:::`._validate_xy`(c("a","b"), c(1,2)),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects non-numeric y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2), c("a","b")),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects unequal lengths", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,2)),
               "x and y must have the same length")
})

test_that("._validate_xy() rejects too-short input", {
  expect_error(mathmodels:::`._validate_xy`(c(1), c(2), min_len = 2),
               "x and y must have length >= 2")
})

test_that("._validate_xy() rejects duplicate x values", {
  expect_error(mathmodels:::`._validate_xy`(c(1,1,2), c(1,2,3)),
               "x must not contain duplicate values")
})

test_that("._validate_xy() rejects NA in x", {
  expect_error(mathmodels:::`._validate_xy`(c(1,NA,3), c(4,5,6)),
               "x and y must not contain NA")
})

test_that("._validate_xy() rejects Inf in y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,Inf,3)),
               "x and y must not contain")
})

# --- interp_linear ---

test_that("interp_linear() returns expected structure", {
  res = interp_linear(c(1,2,3), c(2,4,6), xout = c(1, 1.5, 2, 2.5, 3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "linear")
  expect_equal(nrow(res), 5)
})

test_that("interp_linear() passes through known points", {
  res = interp_linear(c(1,3,5), c(10,30,50), xout = c(1,3,5))
  expect_equal(res$y, c(10, 30, 50))
})

test_that("interp_linear() interpolates midpoints", {
  res = interp_linear(c(0,2), c(0,4), xout = 1)
  expect_equal(res$y, 2)
})

test_that("interp_linear() extrapolates outside range (constant)", {
  res = interp_linear(c(1,2), c(1,2), xout = 0)
  expect_true(is.na(res$y))
})

test_that("interp_linear() validates input via ._validate_xy", {
  expect_error(interp_linear(c(1,1), c(1,2), xout = c(1,2)), "duplicate")
})

test_that("interp_linear() handles single xout", {
  res = interp_linear(c(1,2,3), c(10,20,30), xout = 2)
  expect_equal(nrow(res), 1)
  expect_equal(res$x, 2)
  expect_equal(res$y, 20)
})
