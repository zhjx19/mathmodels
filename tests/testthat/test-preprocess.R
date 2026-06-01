test_that("standardize works", {
  x = c(4, 1, NA, 5, 8)
  r = standardize(x)
  expect_type(r, "double")
  expect_equal(length(r), 5)
  expect_equal(mean(r, na.rm=TRUE), 0, tolerance=1e-8)
  expect_equal(sd(r, na.rm=TRUE), 1, tolerance=1e-8)
})

test_that("normalize works", {
  x = c(3, 4)
  r = normalize(x)
  expect_equal(sqrt(sum(r^2)), 1, tolerance=1e-8)
})

test_that("rescale works", {
  x = c(2, 4, 6, 8)
  r = rescale(x)
  expect_equal(range(r, na.rm=TRUE), c(0, 1))
  r2 = rescale(x, type="-")
  expect_equal(r2, rev(r))
})

test_that("rescale_middle works", {
  PH = 6:9
  r = rescale_middle(PH, 7)
  expect_equal(r[PH==7], 1)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("rescale_interval works", {
  Temp = c(35.2, 35.8, 36.6, 37.1, 37.8, 38.4)
  r = rescale_interval(Temp, 36, 37)
  expect_true(all(r >= 0 & r <= 1, na.rm=TRUE))
  expect_equal(r[Temp >= 36 & Temp <= 37], 1)
})

test_that("rescale_extreme works", {
  x = c(2, 4, 6, 8)
  r = rescale_extreme(x, "+")
  expect_equal(r[4], 1)
  r2 = rescale_extreme(x, "-")
  expect_equal(r2[1], 1)
})

test_that("rescale_initial works", {
  x = c(10, 5, 20)
  r = rescale_initial(x, "+")
  expect_equal(r[1], 1)
  r2 = rescale_initial(x, "-")
  expect_equal(r2[2], 2)
})

test_that("rescale_mean works", {
  x = c(2, 4, 6)
  r = rescale_mean(x)
  expect_equal(mean(r), 1, tolerance=1e-8)
})

test_that("to_positive works", {
  x = c(1, 3, 5)
  r = to_positive(x, "minmax")
  expect_equal(r, c(4, 2, 0))
  r2 = to_positive(x, "reciprocal")
  expect_equal(r2, 1/x)
})

test_that("preprocess input validation", {
  expect_error(standardize("a"), "must be a numeric vector")
  expect_error(rescale("a"), "must be a numeric vector")
  expect_error(rescale(c(1,1,1)), "zero range")
  expect_error(rescale_middle("a", 1), "must be a numeric vector")
})
