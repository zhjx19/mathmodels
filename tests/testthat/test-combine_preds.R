test_that("combine_preds works with 2 values", {
  r = combine_preds(c(100, 102))
  expect_equal(r$w, c(0.5, 0.5))
  expect_equal(r$a, 101)
})

test_that("combine_preds works with 3 values", {
  r = combine_preds(c(100, 102, 98))
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_type(r$a, "double")
})

test_that("combine_preds input validation", {
  expect_error(combine_preds(c(1)), "length >= 2")
})
