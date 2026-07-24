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
