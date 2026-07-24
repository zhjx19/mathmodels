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
