test_that("fuzzy_eval works for all types", {
  w = c(0.3, 0.3, 0.3, 0.1)
  R = matrix(c(0.8, 0.7, 0.6, 0.7,
               0.1, 0.2, 0.2, 0.1,
               0.1, 0.1, 0.2, 0.2), nrow=3, byrow=TRUE)

  for(t in 1:5) {
    r = fuzzy_eval(w, R, t)
    expect_equal(sum(r), 1, tolerance=1e-8)
    expect_length(r, 3)
  }
})

test_that("fuzzy_eval input validation", {
  R2 = matrix(c(0.8, 0.2, 0.3, 0.7, 0.6, 0.4), nrow = 2)
  expect_error(fuzzy_eval("a", R2, 1), "must be a numeric vector")
  expect_error(fuzzy_eval(c(0.3,0.3), R2, 1), "must equal ncol")
  expect_error(fuzzy_eval(c(0.3,0.3,0.4), R2, 6), "must be 1, 2, 3, 4, or 5")
})

test_that("fuzzy_eval rejects NA or out-of-range R", {
  w = c(0.5, 0.5)
  expect_error(fuzzy_eval(w, matrix(c(0.5, NA, 0.3, 0.7), nrow = 2), 3),
               "no NA")
  expect_error(fuzzy_eval(w, matrix(c(0.5, 1.5, 0.3, 0.7), nrow = 2), 3),
               "between 0 and 1")
  expect_error(fuzzy_eval(w, matrix(c(0.5, -0.1, 0.3, 0.7), nrow = 2), 3),
               "between 0 and 1")
})

test_that("fuzzy_eval rejects NA or negative w", {
  R2 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2)
  expect_error(fuzzy_eval(c(0.5, NA), R2, 3), "no NA")
  expect_error(fuzzy_eval(c(-0.5, 1.5), R2, 3), "non-negative")
})
