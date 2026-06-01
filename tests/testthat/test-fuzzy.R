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
  expect_error(fuzzy_eval("a", matrix(1:6, nrow=2), 1), "must be a numeric vector")
  expect_error(fuzzy_eval(c(0.3,0.3), matrix(1:6, nrow=2), 1), "must equal ncol")
  expect_error(fuzzy_eval(c(0.3,0.3,0.4), matrix(1:6, nrow=2), 6), "must be 1, 2, 3, 4, or 5")
})
