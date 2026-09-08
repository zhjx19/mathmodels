# Tests for Rubin pooling (mi_pool.R)

test_that("rubin_pool: m=2 closed-form hand calculation", {
  # Q = (2, 4), U (se^2) = (1, 9)
  # Qbar = 3, Ubar = 5, B = (1+1)/1 = 2, T = 5 + (1 + 1/2) * 2 = 8
  # r = 1.5 * 2 / 5 = 0.6
  # df = (2-1) * (1 + 5/3)^2 = (8/3)^2 = 7.111...
  r = rubin_pool(est = c(2, 4), err = c(1, 3))
  expect_s3_class(r, "rubin_pool")
  expect_equal(r$estimate, 3)
  expect_equal(r$var_within, 5)
  expect_equal(r$var_between, 2)
  expect_equal(r$var_total, 8)
  expect_equal(r$std_error, sqrt(8))
  expect_equal(r$r, 0.6, tolerance = 1e-12)
  expect_equal(r$df, (1 + 5/3)^2, tolerance = 1e-12)
  expect_equal(r$lcl, 3 - stats::qt(0.975, df = (1 + 5/3)^2) * sqrt(8),
               tolerance = 1e-10)
  expect_equal(r$ucl, 3 + stats::qt(0.975, df = (1 + 5/3)^2) * sqrt(8),
               tolerance = 1e-10)
  expect_equal(r$m, 2L)
})

test_that("rubin_pool: fmi / rel_eff follow Rubin definitions", {
  r = rubin_pool(est = c(2, 4), err = c(1, 3))
  # lambda = (r + 2/(df+3)) / (r + 1) with r = 0.6, df = (8/3)^2
  lam = (0.6 + 2 / ((8/3)^2 + 3)) / 1.6
  expect_equal(r$fmi, lam, tolerance = 1e-12)
  expect_equal(r$rel_eff, 1 / (1 + lam / 2), tolerance = 1e-12)
  # t and p-value against theta = 0
  expect_equal(r$t, 3 / sqrt(8), tolerance = 1e-12)
  expect_equal(r$p_value,
               2 * stats::pt(3 / sqrt(8), df = (8/3)^2, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("rubin_pool: identical estimates give B = 0 and infinite df", {
  r = rubin_pool(est = c(1, 1, 1), err = c(2, 2, 2))
  expect_equal(r$var_between, 0)
  expect_equal(r$var_total, r$var_within)
  expect_equal(r$df, Inf)
  # r = 0 -> lambda collapses to r/(r+1) = 0
  expect_equal(r$fmi, 0)
})

test_that("rubin_pool: Barnard-Rubin adjusted df", {
  # edf = 120: gamma = (1 + 1/m) B / T = 1.5 * 2 / 8 = 0.375
  # v_obs = (1 - 0.375) * 120 * 121 / 123
  # df* = 1 / (1/df_rubin + 1/v_obs)
  v_obs = (1 - 0.375) * 120 * 121 / 123
  df_star = 1 / (1 / ((8/3)^2) + 1 / v_obs)
  r = rubin_pool(est = c(2, 4), err = c(1, 3), edf = 120)
  expect_equal(r$df, df_star, tolerance = 1e-12)
  expect_equal(r$edf, 120)
  expect_lt(r$df, v_obs)  # adjusted df never exceeds the v_obs ceiling
})

test_that("rubin_pool: edf = Inf is the default (no adjustment)", {
  r1 = rubin_pool(c(2, 4), c(1, 3))
  r2 = rubin_pool(c(2, 4), c(1, 3), edf = Inf)
  expect_equal(r1$df, r2$df)
})

test_that("rubin_pool: m = 1 degenerates gracefully", {
  r = rubin_pool(est = 2, err = 3)
  expect_equal(r$m, 1L)
  expect_equal(r$var_between, 0)
  expect_equal(r$var_total, 9)
  expect_equal(r$df, Inf)
  expect_equal(r$std_error, 3)
  expect_equal(r$p_value, 2 * stats::pt(2/3, df = Inf, lower.tail = FALSE))
})

test_that("rubin_pool: all-zero se with identical est gives t = Inf", {
  # est all equal AND err all zero -> T = 0, se = 0, delta != 0 -> |t| = Inf
  r = rubin_pool(est = c(2, 2), err = c(0, 0))
  expect_equal(r$std_error, 0)
  expect_equal(r$t, Inf)
  expect_equal(r$p_value, 0)
})

test_that("rubin_pool: input validation", {
  expect_error(rubin_pool("a", 1), "numeric")
  expect_error(rubin_pool(1:3, 1:2), "same length")
  expect_error(rubin_pool(numeric(0), numeric(0)), "at least one")
  expect_error(rubin_pool(c(1, NA), c(1, 1)), "finite")
  expect_error(rubin_pool(c(1, 1), c(-1, 1)), "non-negative")
  expect_error(rubin_pool(1:2, 1:2, alpha = 1.5), "alpha")
  expect_error(rubin_pool(1:2, 1:2, alpha = 0), "alpha")
  expect_error(rubin_pool(1:2, 1:2, edf = -5), "edf")
  expect_error(rubin_pool(1:2, 1:2, edf = 0), "edf")
})
