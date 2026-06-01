test_that("compute_mf works", {
  th = c(0.05, 0.15, 0.25, 0.5)
  r = compute_mf(0.07, th)
  expect_length(r, length(th))
  expect_equal(sum(r), 1, tolerance=1e-8)
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
