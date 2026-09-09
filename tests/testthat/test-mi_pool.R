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
  expect_equal(r$fmi, (0.6 + 2 / (df_star + 3)) / 1.6, tolerance = 1e-12)
  expect_equal(r$rel_eff, 1 / (1 + r$fmi / 2), tolerance = 1e-12)
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

# ---------------------------------------------------------------------------
# pool_fit
# ---------------------------------------------------------------------------

test_that("pool_fit: pools lm coefficients across imputations", {
  set.seed(2)
  d = data.frame(x = rnorm(50), z = rnorm(50))
  d$y = 1 + 2 * d$x + 0.5 * d$z + rnorm(50)
  d$x[c(3, 17)] = NA
  imps = impute_multiple(d, m = 5, maxit = 3, seed = 21)
  fits = lapply(imps, \(dat) lm(y ~ x + z, dat))
  pooled = pool_fit(fits)
  expect_s3_class(pooled, "tbl_df")
  expect_equal(names(pooled),
    c("term", "estimate", "std_error", "lcl", "ucl", "df",
      "t", "p_value", "fmi", "r", "rel_eff"))
  expect_setequal(pooled$term, c("(Intercept)", "x", "z"))
  # cross-check the x row against a manual rubin_pool call
  manual = rubin_pool(
    est = vapply(fits, \(f) unname(coef(f)[2]), numeric(1)),
    err = vapply(fits, \(f) unname(sqrt(diag(vcov(f)))[2]), numeric(1)),
    edf = fits[[1]]$df.residual
  )
  row_x = pooled[pooled$term == "x", ]
  expect_equal(row_x$estimate, manual$estimate, tolerance = 1e-10)
  expect_equal(row_x$std_error, manual$std_error, tolerance = 1e-10)
  expect_equal(row_x$df, manual$df, tolerance = 1e-10)
  # true beta = 2 should fall inside the CI (sanity, not exact)
  expect_lt(row_x$lcl, 2)
  expect_gt(row_x$ucl, 2)
})

test_that("pool_fit: glm works and passes df.residual", {
  set.seed(3)
  n = 80
  x = rnorm(n)
  d = data.frame(x = x)
  d$y = rbinom(n, 1, plogis(-1 + 2 * x))
  d$x[c(5, 33)] = NA
  imps = impute_multiple(d, .cols = x, m = 4, maxit = 3, seed = 8)
  fits = lapply(imps, \(dat) glm(y ~ x, binomial(), dat))
  pooled = pool_fit(fits)
  expect_s3_class(pooled, "tbl_df")
  row_x = pooled[pooled$term == "x", ]
  manual = rubin_pool(
    est = vapply(fits, \(f) unname(coef(f)[2]), numeric(1)),
    err = vapply(fits, \(f) unname(sqrt(diag(vcov(f)))[2]), numeric(1)),
    edf = fits[[1]]$df.residual
  )
  expect_equal(row_x$df, manual$df, tolerance = 1e-10)
})

test_that("pool_fit: NA coefficients are dropped with a warning", {
  d = data.frame(y = 1:6 + rnorm(6), x1 = 1:6, x2 = 2:7)  # collinear
  fits = list(lm(y ~ x1 + x2, d), lm(y ~ x1 + x2, d))
  expect_warning(pooled <- pool_fit(fits), "NA coefficient")
  expect_true(all(c("x1", "(Intercept)") %in% pooled$term))
  expect_false("x2" %in% pooled$term)
})

test_that("pool_fit: mismatched coefficient names error", {
  f1 = lm(mpg ~ wt, mtcars)
  f2 = lm(mpg ~ cyl, mtcars)
  expect_error(pool_fit(list(f1, f2)), "coefficient")
})

test_that("pool_fit: print method runs", {
  set.seed(4)
  d = data.frame(x = rnorm(40), z = rnorm(40))
  d$y = 0.5 * d$x + rnorm(40)
  d$x[7] = NA
  imps = impute_multiple(d, .cols = x, m = 3, maxit = 3, seed = 5)
  fits = lapply(imps, \(dat) lm(y ~ x, dat))
  out = capture.output(print(pool_fit(fits)))
  expect_true(any(grepl("term", out)))
})

test_that("pool_fit: single-coefficient model pools across all fits", {
  # three fits with DIFFERENT estimates: y ~ 0 + x means term = x only
  noise = c(0.02, -0.01, 0.03, -0.02, 0.01)
  d1 = data.frame(x = 1:5, y = 1:5 * 1.0 + noise)
  d2 = data.frame(x = 1:5, y = 1:5 * 2.0 + noise)
  d3 = data.frame(x = 1:5, y = 1:5 * 3.0 + noise)
  fits = list(lm(y ~ 0 + x, d1), lm(y ~ 0 + x, d2), lm(y ~ 0 + x, d3))
  pooled = pool_fit(fits)
  expect_equal(nrow(pooled), 1)
  expect_equal(pooled$term, "x")
  # pooled estimate must be the mean of (~1, ~2, ~3) = ~2, not fit 1's ~1
  expect_equal(pooled$estimate, 2, tolerance = 0.1)
  # and matches manual rubin_pool on the same est/err vectors
  manual = rubin_pool(
    est = vapply(fits, \(f) unname(coef(f)[1]), numeric(1)),
    err = vapply(fits, \(f) unname(sqrt(diag(vcov(f)))[1]), numeric(1)),
    edf = fits[[1]]$df.residual
  )
  expect_equal(pooled$std_error, manual$std_error, tolerance = 1e-10)
  expect_equal(pooled$df, manual$df, tolerance = 1e-10)
})
