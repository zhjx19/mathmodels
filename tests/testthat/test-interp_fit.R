# =============================================================================
# Unit tests for the interpolation & curve fitting module interp_fit.R
# =============================================================================

# --- ._validate_xy ---

test_that("._validate_xy() rejects matrix input", {
  expect_error(mathmodels:::`._validate_xy`(matrix(1:4, 2, 2), c(1, 2, 3, 4)),
               "must be vectors")
})

test_that("._validate_xy() rejects non-numeric x", {
  expect_error(mathmodels:::`._validate_xy`(c("a","b"), c(1,2)),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects non-numeric y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2), c("a","b")),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects unequal lengths", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,2)),
               "x and y must have the same length")
})

test_that("._validate_xy() rejects too-short input", {
  expect_error(mathmodels:::`._validate_xy`(c(1), c(2), min_len = 2),
               "x and y must have length >= 2")
})

test_that("._validate_xy() rejects duplicate x values", {
  expect_error(mathmodels:::`._validate_xy`(c(1,1,2), c(1,2,3)),
               "x must not contain duplicate values")
})

test_that("._validate_xy() rejects NA in x", {
  expect_error(mathmodels:::`._validate_xy`(c(1,NA,3), c(4,5,6)),
               "x and y must not contain NA")
})

test_that("._validate_xy() rejects Inf in y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,Inf,3)),
               "x and y must not contain")
})

# --- interp_linear ---

test_that("interp_linear() returns expected structure", {
  res = interp_linear(c(1,2,3), c(2,4,6), xout = c(1, 1.5, 2, 2.5, 3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "linear")
  expect_equal(nrow(res), 5)
})

test_that("interp_linear() passes through known points", {
  res = interp_linear(c(1,3,5), c(10,30,50), xout = c(1,3,5))
  expect_equal(res$y, c(10, 30, 50))
})

test_that("interp_linear() interpolates midpoints", {
  res = interp_linear(c(0,2), c(0,4), xout = 1)
  expect_equal(res$y, 2)
})

test_that("interp_linear() extrapolates outside range (constant)", {
  res = interp_linear(c(1,2), c(1,2), xout = 0)
  expect_true(is.na(res$y))
})

test_that("interp_linear() validates input via ._validate_xy", {
  expect_error(interp_linear(c(1,1), c(1,2), xout = c(1,2)), "duplicate")
})

# --- interp_spline ---

test_that("interp_spline() returns expected structure", {
  res = interp_spline(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "spline")
})

test_that("interp_spline() passes through known points", {
  res = interp_spline(c(1,3,5,7), c(10,30,50,70), xout = c(1,3,5,7))
  expect_equal(res$y, c(10, 30, 50, 70))
})

test_that("interp_spline() interpolates between points", {
  res = interp_spline(c(0, 1, 2, 3), c(0, 2, 4, 6), xout = 0.5)
  expect_gt(res$y, 0)
  expect_lt(res$y, 4)
})

test_that("interp_spline() with spar uses smooth.spline", {
  res = interp_spline(c(1,2,3,4,5), c(2.1,3.9,5.8,8.2,9.9), xout = c(1,3,5), spar = 0.5)
  expect_identical(attr(res, "method"), "smooth_spline")
  expect_equal(nrow(res), 3)
})

test_that("interp_spline() requires at least 4 points", {
  expect_error(interp_spline(c(1,2,3), c(1,2,3), xout = c(1,2)),
               "length >= 4")
})

test_that("interp_linear() handles single xout", {
  res = interp_linear(c(1,2,3), c(10,20,30), xout = 2)
  expect_equal(nrow(res), 1)
  expect_equal(res$x, 2)
  expect_equal(res$y, 20)
})

# --- interp_poly ---

test_that("interp_poly() returns expected structure", {
  res = interp_poly(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3), degree = 2)
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "poly")
})

test_that("interp_poly() degree=1 fits a straight line", {
  res = interp_poly(c(0,2,4), c(0,4,8), xout = c(0,2,4), degree = 1)
  expect_equal(res$y, c(0, 4, 8))
})

test_that("interp_poly() degree too high for data triggers error", {
  expect_error(interp_poly(c(1,2,3), c(1,2,3), xout = c(1,2,3), degree = 5))
})

test_that("interp_poly() validates min_len based on degree", {
  expect_error(interp_poly(c(1,2), c(1,2), xout = c(1,2), degree = 3),
               "length >= 4")
})

# --- interp_hermite ---

test_that("interp_hermite() returns expected structure", {
  res = interp_hermite(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "hermite")
})

test_that("interp_hermite() passes through known points", {
  res = interp_hermite(c(1,3,5,7), c(10,30,50,70), xout = c(1,3,5,7))
  expect_equal(res$y, c(10, 30, 50, 70))
})

test_that("interp_hermite() preserves monotonicity", {
  x = c(1, 2, 3, 4, 5)
  y = c(1, 2, 4, 7, 11)
  res = interp_hermite(x, y, xout = seq(1, 5, by = 0.1))
  expect_true(all(diff(res$y) >= -1e-10))
})

test_that("interp_hermite() validates input", {
  expect_error(interp_hermite(c(1,1,3), c(1,2,3), xout = c(1,2)), "duplicate")
})

# --- fitting helpers ---

test_that("._compute_fit_stats() returns expected columns", {
  d = data.frame(x = 1:3, y = c(1.1, 2.2, 2.9))
  fit = lm(y ~ x, data = d)
  info = mathmodels:::._compute_fit_stats(d$y, fitted(fit), fit)
  expect_s3_class(info, "tbl_df")
  expect_named(info, c("r.squared", "adj.r.squared", "aic", "bic", "sigma"))
  expect_true(info$r.squared >= 0 && info$r.squared <= 1)
})

test_that("._compute_fit_stats() gives r.squared=1 for perfect fit", {
  d = data.frame(x = 1:3, y = 1:3)
  fit = lm(y ~ x, data = d)
  info = mathmodels:::._compute_fit_stats(d$y, fitted(fit), fit)
  expect_equal(info$r.squared, 1)
  expect_equal(info$sigma, 0)
})

test_that("._tidy_fit_coefs() tidy from lm model", {
  fit = lm(y ~ x, data = data.frame(x = 1:10, y = 2 * (1:10) + rnorm(10, sd = 0.1)))
  coefs = mathmodels:::._tidy_fit_coefs(fit)
  expect_s3_class(coefs, "tbl_df")
  expect_named(coefs, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_equal(coefs$term[1], "(Intercept)")
})

test_that("._build_fit_residuals() returns expected structure", {
  resid_tbl = mathmodels:::._build_fit_residuals(c(1, 2, 3), c(1.1, 2.2, 2.9))
  expect_s3_class(resid_tbl, "tbl_df")
  expect_named(resid_tbl, c(".observed", ".fitted", ".residual"))
  expect_equal(resid_tbl$.observed, c(1, 2, 3))
  expect_equal(resid_tbl$.residual[1], 1 - 1.1)
})

# --- poly_fit ---

test_that("poly_fit() returns expected list structure", {
  set.seed(42)
  x = 1:20
  y = 3 + 2*x + rnorm(20, sd = 1)
  res = poly_fit(x, y, degree = 1)
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_s3_class(res$coefficient, "tbl_df")
  expect_s3_class(res$residuals, "tbl_df")
  expect_equal(res$type, "poly")
})

test_that("poly_fit() degree=2 recovers quadratic coefficients", {
  set.seed(42)
  x = seq(-5, 5, length.out = 30)
  y = 1 + 2*x + 0.5*x^2 + rnorm(30, sd = 0.5)
  res = poly_fit(x, y, degree = 2)
  expect_true(abs(res$coefficient$estimate[3] - 0.5) < 0.2)
})

test_that("poly_fit() validates input", {
  expect_error(poly_fit(c(1,2), c(1,2), degree = 3), "length >= 4")
})

test_that("poly_fit() residuals match data length", {
  x = 1:10
  y = 2 + 3*x + rnorm(10, sd = 0.5)
  res = poly_fit(x, y, degree = 1)
  expect_equal(nrow(res$residuals), 10)
  expect_equal(nrow(res$input), 10)
})

test_that("poly_fit() model_info AIC matches stats::AIC", {
  x = 1:20
  y = 3 + 2*x + rnorm(20, sd = 1)
  res = poly_fit(x, y, degree = 1)
  expect_equal(res$model_info$aic, round(stats::AIC(res$model), 4))
  expect_equal(res$model_info$bic, round(stats::BIC(res$model), 4))
})

# --- curve_fit ---

test_that("curve_fit() returns expected list structure", {
  set.seed(42)
  x = 1:20
  y = 2 * exp(0.15 * x) * exp(rnorm(20, sd = 0.1))
  res = curve_fit(x, y, type = "exp")
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_equal(res$type, "exp")
})

test_that("curve_fit(type='exp') recovers parameters", {
  set.seed(42)
  x = seq(0, 5, length.out = 50)
  y = 2 * exp(0.5 * x) * exp(rnorm(50, sd = 0.05))
  res = curve_fit(x, y, type = "exp")
  # coefficient on transformed scale: intercept ~ log(a), slope ~ b
  expect_true(abs(res$coefficient$estimate[1] - log(2)) < 0.3)
  expect_true(abs(res$coefficient$estimate[2] - 0.5) < 0.1)
})

test_that("curve_fit(type='power') recovers parameters", {
  set.seed(42)
  x = seq(1, 10, length.out = 50)
  y = 3 * x^2 * exp(rnorm(50, sd = 0.1))
  res = curve_fit(x, y, type = "power")
  expect_true(abs(exp(res$coefficient$estimate[1]) - 3) < 1)
  expect_true(abs(res$coefficient$estimate[2] - 2) < 0.3)
})

test_that("curve_fit(type='log') recovers parameters", {
  set.seed(42)
  x = seq(1, 20, length.out = 50)
  y = 1 + 3 * log(x) + rnorm(50, sd = 0.3)
  res = curve_fit(x, y, type = "log")
  expect_true(abs(res$coefficient$estimate[1] - 1) < 0.5)
  expect_true(abs(res$coefficient$estimate[2] - 3) < 0.3)
})

test_that("curve_fit(type='hyperbolic') recovers parameters", {
  set.seed(42)
  x = seq(1, 20, length.out = 50)
  y = 5 + 10/x + rnorm(50, sd = 0.3)
  res = curve_fit(x, y, type = "hyperbolic")
  expect_true(abs(res$coefficient$estimate[1] - 5) < 1)
  expect_true(abs(res$coefficient$estimate[2] - 10) < 2)
})

test_that("curve_fit() basic validation for exp type", {
  # y must be > 0 for exp fit
  expect_error(curve_fit(1:5, c(-1, 1, 2, 3, 4), type = "exp"), "positive")
})

test_that("curve_fit() basic validation for power type", {
  # x must be > 0 for power fit
  expect_error(curve_fit(c(1, -2, 3), c(1, 2, 3), type = "power"), "positive")
})

test_that("curve_fit() rejects invalid type", {
  expect_error(curve_fit(1:5, 1:5, type = "unknown"), "should be one of")
})

# --- growth_fit ---

test_that("growth_fit() returns expected list structure", {
  set.seed(42)
  x = seq(0, 10, length.out = 30)
  y = 100 / (1 + exp(-1.5 * (x - 5))) + rnorm(30, sd = 2)
  res = growth_fit(x, y, type = "logistic")
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_equal(res$type, "logistic")
  expect_equal(nrow(res$residuals), 30)
})

test_that("growth_fit(type='logistic') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  L_true = 100; k_true = 1.5; x0_true = 5
  y = L_true / (1 + exp(-k_true * (x - x0_true))) + rnorm(50, sd = 2)
  res = growth_fit(x, y, type = "logistic")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - L_true) < 10)
  expect_true(abs(est[2] - k_true) < 0.5)
  expect_true(abs(est[3] - x0_true) < 1)
})

test_that("growth_fit(type='gompertz') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  L_true = 100; k_true = 0.5; x0_true = 4
  y = L_true * exp(-exp(-k_true * (x - x0_true))) + rnorm(50, sd = 2)
  res = growth_fit(x, y, type = "gompertz")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - L_true) < 10)
  expect_true(abs(est[2] - k_true) < 0.3)
  expect_true(abs(est[3] - x0_true) < 1.5)
})

test_that("growth_fit(type='saturation') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  a_true = 50; b_true = 0.5
  y = a_true * (1 - exp(-b_true * x)) + rnorm(50, sd = 1)
  res = growth_fit(x, y, type = "saturation")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - a_true) < 5)
  expect_true(abs(est[2] - b_true) < 0.2)
})

test_that("growth_fit(type='mm') recovers parameters", {
  set.seed(42)
  x = seq(0.5, 20, length.out = 50)
  a_true = 30; b_true = 3
  y = a_true * x / (b_true + x) + rnorm(50, sd = 0.5)
  res = growth_fit(x, y, type = "mm")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - a_true) < 5)
  expect_true(abs(est[2] - b_true) < 2)
})

test_that("growth_fit() rejects invalid type", {
  expect_error(growth_fit(1:10, 1:10, type = "unknown"), "should be one of")
})

test_that("growth_fit() validates input", {
  expect_error(growth_fit(c(1,2), c(1,2), type = "logistic"), "length >= 3")
})
