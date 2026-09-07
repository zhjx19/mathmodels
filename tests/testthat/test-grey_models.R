test_that("GM11 works", {
  x = c(100, 120, 145, 175, 210)
  r = GM11(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u", "lambda", "rng"))
  expect_length(r$fitted, length(x))
  expect_type(r$pnext, "double")
  expect_type(r$f, "closure")
})

test_that("GM11 input validation", {
  expect_error(GM11("a"), "must be a numeric vector")
  expect_error(GM11(c(1,2,3)), "length >= 4")
})

test_that("GM11 rejects non-positive, missing, and constant series", {
  expect_error(GM11(c(-1, 2, 3, 4)), "strictly positive")
  expect_error(GM11(c(0, 2, 3, 4)), "strictly positive")
  expect_error(GM11(c(1, NA, 3, 4)), "strictly positive")
  expect_error(GM11(c(5, 5, 5, 5)), "essentially constant")
})

test_that("GM11 reproduces the textbook GM(1,1) least-squares fit", {
  x0 = c(15.2, 16.1, 17.9, 18.7, 19.6, 21.3)
  r = GM11(x0)

  # Reference fit: standard GM(1,1) least-squares system solved inline
  n = length(x0)
  x1 = cumsum(x0)
  z = (x1[1:(n - 1)] + x1[2:n]) / 2
  ref = stats::lm.fit(cbind(-z, rep(1, n - 1)), x0[2:n])
  a = as.numeric(ref$coefficients[1])
  b = as.numeric(ref$coefficients[2])

  x1h = (x0[1] - b / a) * exp(-a * (1:n - 1)) + b / a
  fitted_ref = c(x1h[1], diff(x1h))
  pnext_ref = (x0[1] - b / a) * (1 - exp(a)) * exp(-a * n)

  expect_equal(r$fitted[1], x0[1])
  expect_equal(as.numeric(r$u), c(a, b), tolerance = 1e-6)
  expect_equal(r$fitted, fitted_ref, tolerance = 1e-6)
  expect_equal(r$pnext, pnext_ref, tolerance = 1e-6)
})

test_that("DGM21 works", {
  x = c(2.874, 3.278, 3.39, 3.679, 3.77, 3.8)
  r = DGM21(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u"))
  expect_length(r$fitted, length(x))
})

test_that("verhulst works", {
  x = c(4.93, 2.33, 3.87, 4.35, 6.63, 7.15, 5.37, 6.39, 7.81, 8.35)
  r = verhulst(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u"))
  expect_length(r$fitted, length(x))
})

test_that("GM1N works", {
  data = data.frame(
    factor1 = c(50, 55, 60, 65, 70),
    factor2 = c(20, 22, 25, 28, 30),
    output = c(100, 120, 145, 175, 210)
  )
  r = GM1N(data)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pred", "C", "P", "f"))
  expect_length(r$fitted, nrow(data))
  expect_null(r$pred)
  expect_type(r$C, "double")
  expect_type(r$P, "double")
  expect_true(r$C > 0)
  expect_true(r$P >= 0 && r$P <= 1)
})

test_that("GM1N input validation", {
  expect_error(GM1N("a"), "must be a data frame or matrix")
  bad_data = data.frame(factor1 = c(1, NA, 3, 4), output = c(2, 3, 4, 5))
  expect_error(GM1N(bad_data), "strictly positive")
  bad_data2 = data.frame(factor1 = c(1, -2, 3, 4), output = c(2, 3, 4, 5))
  expect_error(GM1N(bad_data2), "strictly positive")
})
