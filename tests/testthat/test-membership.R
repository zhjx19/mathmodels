test_that("tri_mf works", {
  x = 0:10
  r = tri_mf(x, c(3, 6, 8))
  expect_length(r, 11)
  expect_equal(r[7], 1, tolerance=1e-8)
  expect_equal(r[1], 0)
  expect_equal(r[11], 0)
})

test_that("trap_mf works", {
  x = 0:10
  r = trap_mf(x, c(1, 5, 7, 8))
  expect_equal(r[6], 1)
  expect_equal(r[1], 0)
})

test_that("gauss_mf works", {
  r = gauss_mf(0:10, c(2, 5))
  expect_equal(r[6], 1, tolerance=1e-8)
})

test_that("gbell_mf works", {
  r = gbell_mf(0:10, c(2, 4, 6))
  expect_true(all(r >= 0 & r <= 1))
})

test_that("gauss2mf works", {
  r = gauss2mf(0:10, c(1, 3, 3, 4))
  expect_true(all(r >= 0 & r <= 1))
})

test_that("sigmoid_mf works", {
  r = sigmoid_mf(-5:5, c(2, 0))
  expect_equal(r[6], 0.5, tolerance=1e-8)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("dsigmoid_mf works", {
  r = dsigmoid_mf(0:10, c(5, 2, 5, 7))
  expect_true(all(r <= 1))
})

test_that("psigmoid_mf works", {
  r = psigmoid_mf(0:10, c(2, 3, -5, 8))
  expect_true(all(r >= 0 & r <= 1, na.rm=TRUE))
})

test_that("z_mf works", {
  r = z_mf(0:10, c(3, 7))
  expect_equal(r[1], 1)
  expect_equal(r[11], 0)
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("z_mf two-stage spline correctness", {
  # At midpoint (a+b)/2 = 5, value should be exactly 0.5
  r = z_mf(5, c(3, 7))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Monotonic decreasing
  x = seq(3, 7, length.out = 21)
  r = z_mf(x, c(3, 7))
  expect_true(all(diff(r) <= 0))
})

test_that("pi_mf works", {
  r = pi_mf(0:10, c(1, 4, 5, 10))
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("pi_mf two-stage spline correctness", {
  a = 1; b = 4; c = 5; d = 10
  # Midpoint of rising edge: (a+b)/2 = 2.5 → μ = 0.5
  r = pi_mf(2.5, c(a, b, c, d))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Midpoint of falling edge: (c+d)/2 = 7.5 → μ = 0.5
  r = pi_mf(7.5, c(a, b, c, d))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Plateau
  expect_equal(pi_mf(b, c(a,b,c,d)), 1)
  expect_equal(pi_mf(c, c(a,b,c,d)), 1)
  # Monotonic rising
  x = seq(a, b, length.out = 21)
  r = pi_mf(x, c(a, b, c, d))
  expect_true(all(diff(r) >= 0))
  # Monotonic falling
  x = seq(c, d, length.out = 21)
  r = pi_mf(x, c(a, b, c, d))
  expect_true(all(diff(r) <= 0))
})

test_that("s_mf works", {
  r = s_mf(0:10, c(1, 8))
  expect_equal(r[1], 0)
  expect_equal(r[11], 1)
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("s_mf two-stage spline correctness", {
  # At midpoint (a+b)/2 = 4.5, value should be exactly 0.5
  r = s_mf(4.5, c(1, 8))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Monotonic increasing
  x = seq(1, 8, length.out = 21)
  r = s_mf(x, c(1, 8))
  expect_true(all(diff(r) >= 0))
})

test_that("membership input validation", {
  expect_error(tri_mf("a", c(3,6,8)), "must be a numeric vector")
  expect_error(tri_mf(0:10, c(3,6)), "length 3")
})
