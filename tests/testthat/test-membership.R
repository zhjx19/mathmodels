test_that("mf_tri works", {
  x = 0:10
  r = mf_tri(x, c(3, 6, 8))
  expect_length(r, 11)
  expect_equal(r[7], 1, tolerance=1e-8)
  expect_equal(r[1], 0)
  expect_equal(r[11], 0)
})

test_that("mf_trap works", {
  x = 0:10
  r = mf_trap(x, c(1, 5, 7, 8))
  expect_equal(r[6], 1)
  expect_equal(r[1], 0)
})

test_that("mf_gauss works", {
  r = mf_gauss(0:10, c(2, 5))
  expect_equal(r[6], 1, tolerance=1e-8)
})

test_that("mf_gbell works", {
  r = mf_gbell(0:10, c(2, 4, 6))
  expect_true(all(r >= 0 & r <= 1))
})

test_that("mf_gauss2 works", {
  r = mf_gauss2(0:10, c(1, 3, 3, 4))
  expect_true(all(r >= 0 & r <= 1))
})

test_that("mf_sigmoid works", {
  r = mf_sigmoid(-5:5, c(2, 0))
  expect_equal(r[6], 0.5, tolerance=1e-8)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("mf_dsigmoid works", {
  r = mf_dsigmoid(0:10, c(5, 2, 5, 7))
  expect_true(all(r <= 1))
})

test_that("mf_psigmoid works", {
  r = mf_psigmoid(0:10, c(2, 3, -5, 8))
  expect_true(all(r >= 0 & r <= 1, na.rm=TRUE))
})

test_that("mf_z works", {
  r = mf_z(0:10, c(3, 7))
  expect_equal(r[1], 1)
  expect_equal(r[11], 0)
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("mf_z two-stage spline correctness", {
  # At midpoint (a+b)/2 = 5, value should be exactly 0.5
  r = mf_z(5, c(3, 7))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Monotonic decreasing
  x = seq(3, 7, length.out = 21)
  r = mf_z(x, c(3, 7))
  expect_true(all(diff(r) <= 0))
})

test_that("mf_pi works", {
  r = mf_pi(0:10, c(1, 4, 5, 10))
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("mf_pi two-stage spline correctness", {
  a = 1; b = 4; c = 5; d = 10
  # Midpoint of rising edge: (a+b)/2 = 2.5 → μ = 0.5
  r = mf_pi(2.5, c(a, b, c, d))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Midpoint of falling edge: (c+d)/2 = 7.5 → μ = 0.5
  r = mf_pi(7.5, c(a, b, c, d))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Plateau
  expect_equal(mf_pi(b, c(a,b,c,d)), 1)
  expect_equal(mf_pi(c, c(a,b,c,d)), 1)
  # Monotonic rising
  x = seq(a, b, length.out = 21)
  r = mf_pi(x, c(a, b, c, d))
  expect_true(all(diff(r) >= 0))
  # Monotonic falling
  x = seq(c, d, length.out = 21)
  r = mf_pi(x, c(a, b, c, d))
  expect_true(all(diff(r) <= 0))
})

test_that("mf_s works", {
  r = mf_s(0:10, c(1, 8))
  expect_equal(r[1], 0)
  expect_equal(r[11], 1)
  expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("mf_s two-stage spline correctness", {
  # At midpoint (a+b)/2 = 4.5, value should be exactly 0.5
  r = mf_s(4.5, c(1, 8))
  expect_equal(r, 0.5, tolerance = 1e-8)
  # Monotonic increasing
  x = seq(1, 8, length.out = 21)
  r = mf_s(x, c(1, 8))
  expect_true(all(diff(r) >= 0))
})

test_that("membership input validation", {
  expect_error(mf_tri("a", c(3,6,8)), "must be a numeric vector")
  expect_error(mf_tri(0:10, c(3,6)), "length 3")
})

test_that("membership degenerate parameters are rejected", {
  expect_error(mf_tri(0:10, c(2, 2, 8)), "a < b < c")
  expect_error(mf_tri(0:10, c(2, 8, 8)), "a < b < c")
  expect_error(mf_trap(0:10, c(2, 2, 6, 8)), "a < b <= c < d")
  expect_error(mf_trap(0:10, c(1, 3, 7, 7)), "a < b <= c < d")
  expect_error(mf_gauss(0:10, c(0, 5)), "sigma must be positive")
  expect_error(mf_gauss(0:10, c(-1, 5)), "sigma must be positive")
  expect_error(mf_gbell(0:10, c(0, 2, 5)), "a must be positive")
  expect_error(mf_z(0:10, c(3, 3)), "a < b")
  expect_error(mf_s(0:10, c(3, 3)), "a < b")
  expect_error(mf_pi(0:10, c(4, 4, 5, 10)), "a < b <= c < d")
})

test_that("mf_trap one-sided Inf shoulders still work", {
  r1 = mf_trap(0:10, c(-Inf, -Inf, 5, 8))
  expect_equal(r1[1], 1)
  expect_equal(r1[5], 1)
  expect_equal(r1[11], 0)
  r2 = mf_trap(0:10, c(1, 3, Inf, Inf))
  expect_equal(r2[1], 0)
  expect_equal(r2[4], 1)
  expect_equal(r2[11], 1)
})
