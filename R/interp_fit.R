# =============================================================================
# Interpolation and curve fitting for mathmodels
# Implements: linear/poly/spline/hermite interpolation,
#             polynomial fit, linearizable curve fit, nonlinear growth fit
# =============================================================================

# -----------------------------------------------------------------------------
# Internal Helpers
# -----------------------------------------------------------------------------

# Validate x, y input pair
._validate_xy = function(x, y, min_len = 2) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.", call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  if (length(x) < min_len) {
    stop(sprintf("x and y must have length >= %d.", min_len), call. = FALSE)
  }
  if (anyDuplicated(x) != 0) {
    stop("x must not contain duplicate values.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must not contain NA, Inf, or NaN.", call. = FALSE)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Fitting helpers
# -----------------------------------------------------------------------------

# Compute fit stats on the original scale
._compute_fit_stats = function(observed, fitted, df_residual) {
  ss_res = sum((observed - fitted)^2)
  ss_tot = sum((observed - mean(observed))^2)
  r_sq   = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  n      = length(observed)
  p      = n - df_residual
  adj_r_sq = if (ss_tot > 0 && df_residual > 0) {
    1 - (1 - r_sq) * (n - 1) / df_residual
  } else NA_real_
  sigma = if (df_residual > 0) sqrt(ss_res / df_residual) else 0
  nll   = -sum(stats::dnorm(observed, mean = fitted, sd = sigma, log = TRUE))
  aic   = 2 * p + 2 * nll
  bic   = log(n) * p + 2 * nll
  tibble::tibble(
    r.squared     = round(r_sq, 4),
    adj.r.squared = round(adj_r_sq, 4),
    aic           = round(aic, 4),
    bic           = round(bic, 4),
    sigma         = round(sigma, 4)
  )
}

# Build tidy coefficient tibble from lm or nls model
._tidy_fit_coefs = function(model, level = 0.95) {
  coefs = stats::coef(model)
  se    = sqrt(diag(stats::vcov(model)))
  df_residual = if (inherits(model, "lm")) {
    model$df.residual
  } else {
    length(stats::residuals(model)) - length(coefs)
  }
  t_crit = stats::qt((1 - level) / 2, df_residual, lower.tail = FALSE)
  tibble::tibble(
    term      = names(coefs),
    estimate  = as.numeric(coefs),
    std.error = as.numeric(se),
    statistic = as.numeric(coefs) / as.numeric(se),
    p.value   = 2 * stats::pt(abs(as.numeric(coefs) / as.numeric(se)),
                               df_residual, lower.tail = FALSE),
    conf.low  = as.numeric(coefs) - t_crit * as.numeric(se),
    conf.high = as.numeric(coefs) + t_crit * as.numeric(se)
  )
}

# Build residuals tibble
._build_fit_residuals = function(observed, fitted) {
  tibble::tibble(
    .observed = observed,
    .fitted   = fitted,
    .residual = observed - fitted
  )
}

# Assemble unified return value
._assemble_return = function(model, observed, fitted, formula, type, x, y) {
  coefs  = ._tidy_fit_coefs(model)
  df_residual = if (inherits(model, "lm")) {
    model$df.residual
  } else {
    length(fitted) - length(stats::coef(model))
  }
  info   = ._compute_fit_stats(observed, fitted, df_residual)
  resid  = ._build_fit_residuals(observed, fitted)
  list(
    model       = model,
    coefficient = coefs,
    model_info  = info,
    residuals   = resid,
    formula     = formula,
    type        = type,
    input       = tibble::tibble(x = x, y = y)
  )
}

# -----------------------------------------------------------------------------
# Interpolation
# -----------------------------------------------------------------------------

#' Linear Interpolation
#'
#' Performs linear (piecewise) interpolation at specified points.
#'
#' @param x Numeric vector of observed x-values (no duplicates).
#' @param y Numeric vector of observed y-values, same length as `x`.
#' @param xout Numeric vector of x-values where interpolation is desired.
#'
#' @return A tibble with columns `x` and `y`, containing the interpolated
#'   values at `xout`. The attribute `method` is set to `"linear"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_linear(x, y, xout = c(1.5, 2.5, 3.5))
#'
#' @export
interp_linear = function(x, y, xout) {
  ._validate_xy(x, y, min_len = 2)
  yout = stats::approx(x, y, xout = xout)$y
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = "linear"
  res
}

#' Cubic Spline Interpolation
#'
#' Performs cubic spline interpolation. When `spar` is `NULL`, uses
#' `stats::spline()` with the FMM method. Otherwise uses
#' `stats::smooth.spline()` for smoothing.
#'
#' @param x Numeric vector of observed x-values (no duplicates, length >= 4).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#' @param spar Smoothing parameter passed to `smooth.spline()`. `NULL` for
#'   standard (non-smoothing) spline.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set to
#'   `"spline"` or `"smooth_spline"` depending on the variant used.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_spline(x, y, xout = c(1.5, 2.5, 3.5))
#' interp_spline(x, y, xout = c(1.5, 2.5, 3.5), spar = 0.8)
#'
#' @export
interp_spline = function(x, y, xout, spar = NULL) {
  ._validate_xy(x, y, min_len = 4)
  if (is.null(spar)) {
    yout = stats::spline(x, y, xout = xout, method = "fmm")$y
    method = "spline"
  } else {
    fit = stats::smooth.spline(x, y, spar = spar)
    yout = as.numeric(predict(fit, xout)$y)
    method = "smooth_spline"
  }
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = method
  res
}

#' Polynomial Interpolation
#'
#' Fits a polynomial of specified degree via `lm()` and evaluates it at
#' `xout`.
#'
#' @param x Numeric vector of observed x-values (no duplicates).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#' @param degree Integer, degree of the polynomial. Minimum 1.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set
#'   to `"poly"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_poly(x, y, xout = c(1.5, 2.5, 3.5), degree = 3)
#'
#' @export
interp_poly = function(x, y, xout, degree = 3) {
  ._validate_xy(x, y, min_len = degree + 1)
  fit = lm(y ~ poly(x, degree, raw = TRUE))
  yout = predict(fit, newdata = data.frame(x = xout))
  res = tibble::tibble(x = xout, y = as.numeric(yout))
  attr(res, "method") = "poly"
  res
}

#' Piecewise Cubic Hermite Interpolation
#'
#' Shape-preserving (monotone) piecewise cubic Hermite interpolation using
#' the Fritsch-Carlson method via `splinefun(method = "monoH.FC")`.
#'
#' @param x Numeric vector of observed x-values (no duplicates, length >= 2).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set
#'   to `"hermite"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(1, 2, 4, 7, 11)
#' interp_hermite(x, y, xout = c(1.5, 2.5, 3.5))
#'
#' @export
interp_hermite = function(x, y, xout) {
  ._validate_xy(x, y, min_len = 2)
  sf = stats::splinefun(x, y, method = "monoH.FC")
  yout = sf(xout)
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = "hermite"
  res
}

# -----------------------------------------------------------------------------
# Curve Fitting
# -----------------------------------------------------------------------------

#' Polynomial Fit
#'
#' Fits a polynomial regression model via `lm()`.
#'
#' @param x Numeric vector of observed predictor values.
#' @param y Numeric vector of observed response values.
#' @param degree Integer, degree of the polynomial. Default is 1 (linear).
#'
#' @return A named list with elements: `model` (lm object), `coefficient`
#'   (tibble with estimates, SE, t-value, p-value, 95% CI), `model_info`
#'   (tibble with r.squared, adj.r.squared, AIC, BIC, sigma), `residuals`
#'   (tibble of .observed, .fitted, .residual), `formula` (character),
#'   `type` (character), and `input` (tibble of x, y).
#'
#' @examples
#' x = 1:20
#' y = 3 + 2*x + rnorm(20, sd = 2)
#' poly_fit(x, y, degree = 1)
#' poly_fit(x, y, degree = 2)
#'
#' @export
poly_fit = function(x, y, degree = 1) {
  ._validate_xy(x, y, min_len = degree + 1)
  fit = lm(y ~ poly(x, degree, raw = TRUE))
  fitted_vals = as.numeric(stats::fitted(fit))
  formula_str = paste0("y ~ poly(x, ", degree, ", raw = TRUE)")
  ._assemble_return(fit, y, fitted_vals, formula_str, "poly", x, y)
}
