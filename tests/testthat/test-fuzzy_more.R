test_that("compute_mf tri default is backward-compatible", {
  th = c(0.05, 0.15, 0.25, 0.5)
  r = compute_mf(0.07, th)
  expect_length(r, length(th))
  expect_equal(r, c(0.8, 0.2, 0, 0), tolerance = 1e-8)
})

test_that("compute_mf tri builder explicit", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0.07, th, .builder = "tri")
  expect_length(r, 3)
  expect_equal(sum(r), 1, tolerance = 1e-8)
})

test_that("compute_mf edge case: below first threshold", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0.01, th)
  expect_equal(r[1], 1)
  expect_equal(sum(r[-1]), 0)
})

test_that("compute_mf_funs works", {
  th = c(0.05, 0.15, 0.25)
  mfs = compute_mf_funs(th)
  expect_length(mfs, length(th))
  expect_type(mfs[[1]], "closure")
})

test_that("compute_mf gauss builder default sigma", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0, th, .builder = "gauss")
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(r[1], 1)  # far left → fully level 1
})

test_that("compute_mf gauss builder custom sigma", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0.15, th, .builder = "gauss", sigma = 0.05)
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
  # At center of level 2, membership peaks there
  expect_equal(r[2], 1, tolerance = 1e-6)
})

test_that("compute_mf gauss builder far right", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(1, th, .builder = "gauss", sigma = 0.1)
  expect_length(r, 3)
  expect_equal(r[3], 1, tolerance = 1e-6)
})

test_that("compute_mf sigmoid builder default slope", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0, th, .builder = "sigmoid")
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("compute_mf sigmoid builder custom slope", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(0.15, th, .builder = "sigmoid", slope = 50)
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("compute_mf sigmoid builder far right", {
  th = c(0.05, 0.15, 0.25)
  r = compute_mf(1, th, .builder = "sigmoid", slope = 50)
  expect_true(r[3] > 0.99)
})

test_that("compute_mf custom builder function", {
  th = c(0.05, 0.15, 0.25)
  exp_decay = function(knots, rate = 10) {
    n = length(knots)
    lapply(seq_len(n), function(i) {
      force(i)
      function(x) exp(-rate * abs(x - knots[i]))
    })
  }
  r = compute_mf(0.15, th, .builder = exp_decay, rate = 10)
  expect_length(r, 3)
  expect_true(all(r >= 0 & r <= 1))
  expect_equal(r[2], 1, tolerance = 1e-6)
})

test_that("compute_mf_funs custom builder returns closures", {
  th = c(0.05, 0.15, 0.25)
  custom = function(knots, ...) {
    n = length(knots)
    lapply(seq_len(n), function(i) {
      force(i)
      function(x) 1 / (1 + abs(x - knots[i]))
    })
  }
  mfs = compute_mf_funs(th, .builder = custom)
  expect_length(mfs, 3)
  expect_type(mfs[[1]], "closure")
})

test_that("compute_mf input validation for .builder", {
  th = c(0.05, 0.15, 0.25)
  expect_error(compute_mf(0.1, th, .builder = "unknown"),
               "should be one of")
  expect_error(compute_mf(0.1, th, .builder = 1L),
               "must be a character string or a function")
})

test_that("compute_mf_funs .builder = gauss respects sigma validation", {
  expect_error(compute_mf_funs(c(0, 1), .builder = "gauss", sigma = -1),
               "sigma must be a positive")
})

test_that("defuzzify weighted_average", {
  mu = c(0.318, 0.351, 0.203, 0.128)
  scores = c(30, 60, 75, 90)
  r = defuzzify(mu, scores, "weighted_average")
  expect_type(r, "double")
  expect_equal(r, sum(mu*scores))
})

test_that("defuzzify max_membership", {
  r = defuzzify(c(0.1, 0.3, 0.5, 0.1), c(10, 20, 30, 40), "max_membership")
  expect_equal(r, 30)
})

test_that("defuzzify centroid", {
  mu = c(0.318, 0.351, 0.203, 0.128)
  scores = c(30, 60, 75, 90)
  r = defuzzify(mu, scores, "centroid")
  expect_equal(r, sum(mu*scores)/sum(mu))
})

test_that("defuzzify input validation", {
  expect_error(defuzzify("a", c(1,2,3)), "must be a numeric vector")
  expect_error(defuzzify(c(1,2), c(3,4,5)), "same length")
})

test_that("defuzzify unknown method errors", {
  expect_error(defuzzify(c(0.1, 0.3, 0.5, 0.1), c(10, 20, 30, 40), "bogus"),
               "method must be one of")
})

test_that("defuzzify max_membership ties return the first maximum", {
  r = defuzzify(c(0.5, 0.5, 0.0, 0.1), c(10, 20, 30, 40), "max_membership")
  expect_length(r, 1)
  expect_equal(r, 10)
})

test_that("fuzzy_eval type 5 tolerates a grade with zero membership", {
  w = c(0.2, 0.1, 0.15, 0.3, 0.25)
  R = matrix(c(0,   0,   0,   0,   0,     # grade 1: no factor supports it
               0.3, 0.4, 0.5, 0.2, 0.6,
               0.7, 0.6, 0.5, 0.8, 0.4), nrow = 3, byrow = TRUE)

  r = fuzzy_eval(w, R, 5)

  expect_length(r, 3)
  expect_false(anyNA(r))
  expect_equal(sum(r), 1, tolerance = 1e-12)
  expect_equal(r[1], 0)          # impossible grade stays 0
  # regression: this used to raise "missing value where TRUE/FALSE needed"
  # from 0/0 = NaN, never reaching the function's own error message
  expect_no_error(fuzzy_eval(w, R, 5))
})

test_that("fuzzy_eval errors cleanly when every grade is impossible", {
  R0 = matrix(0, nrow = 2, ncol = 2)
  expect_error(fuzzy_eval(c(0.5, 0.5), R0, 5),
               "Normalization failed: sum of B is zero")
})

test_that("fuzzy_eval rejects negative R and points at mf_dsigmoid", {
  expect_error(fuzzy_eval(c(0.5, 0.5), matrix(c(-0.2, 0.5, 0.5, 1), 2, 2), 3),
               "mf_dsigmoid")
})

test_that("mf_dsigmoid range depends on params (documented as [-1, 1])", {
  x = seq(0, 10, by = 0.01)
  expect_true(all(mf_dsigmoid(x, c(5, 2, 5, 7)) >= -1e-9))  # c1 < c2 -> non-negative
  expect_true(any(mf_dsigmoid(x, c(5, 7, 5, 2)) < 0))       # c1 > c2 -> negative
  expect_true(all(abs(mf_dsigmoid(x, c(5, 7, 5, 2))) <= 1 + 1e-9))
})
