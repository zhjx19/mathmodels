# ============================================================================
# Shared internal input checks for the evaluation family.
#
# A predicate plus an rlang error pipe: the *message* is supplied verbatim by
# the caller, so every pre-existing error string is preserved byte for byte.
# Only the error class and the reported call site change (base `error` ->
# `rlang_error`, attributed to the calling exported function via `call`).
#
# Naming follows the repo-wide convention for internal helpers (dot prefix);
# see `.mv_matrix()` / `.check_numeric_cols()` for the same idea in other
# modules.
# ============================================================================

#' Check a matrix- or data-frame-like argument
#'
#' @param x Object to check.
#' @param msg Error message, passed through verbatim.
#' @param ... Forwarded to [rlang::abort()].
#' @param df_only If `TRUE`, only a data frame is accepted.
#' @param coerce If `TRUE`, convert an accepted matrix to a data frame.
#' @param arg,call Argument name and caller environment for error reporting.
#' @return `x`, or its `as.data.frame()` when `coerce = TRUE`.
#' @noRd
.check_matrix_like = function(x, msg, ..., df_only = FALSE, coerce = FALSE,
                              arg = rlang::caller_arg(x),
                              call = rlang::caller_env()) {
  ok = if (df_only) is.data.frame(x) else (is.data.frame(x) || is.matrix(x))
  if (!ok) rlang::abort(msg, ..., arg = arg, call = call)
  if (coerce && is.matrix(x)) return(as.data.frame(x))
  x
}

#' Check and canonicalise indicator directions
#'
#' @param index Character vector, or `NULL` for an all-default vector.
#' @param m Required length (number of indicators).
#' @param msg_length Error message for a length mismatch.
#' @param msg_allowed Error message for a disallowed value; `NULL` skips the check.
#' @param ... Forwarded to [rlang::abort()].
#' @param allow_na Whether `NA` is an allowed direction.
#' @param default Value used when `index` is `NULL`.
#' @param arg,call Argument name and caller environment for error reporting.
#' @return A character vector of length `m`.
#' @noRd
.check_index = function(index, m, msg_length, msg_allowed = NULL, ...,
                        allow_na = TRUE, default = NA_character_,
                        arg = rlang::caller_arg(index),
                        call = rlang::caller_env()) {
  if (is.null(index)) index = rep(default, m)
  if (length(index) != m) rlang::abort(msg_length, ..., arg = arg, call = call)
  if (!is.null(msg_allowed)) {
    ok = if (allow_na) {
      all(is.na(index) | index %in% c("+", "-"))
    } else {
      all(index %in% c("+", "-"))
    }
    if (!ok) rlang::abort(msg_allowed, ..., arg = arg, call = call)
  }
  index
}

#' Rescale indicators by direction, for data frames and matrices alike
#'
#' @param X Data frame or matrix of indicators (one per column).
#' @param index Character vector of `"+"`, `"-"` or `NA`, length `ncol(X)`.
#' @param a,b Target range passed through to [rescale()].
#' @return A data frame with the same dimensions.
#' @noRd
.rescale_indicators = function(X, index, a = 0, b = 1) {
  X = as.data.frame(X)
  for (j in which(index == "+")) X[[j]] = rescale(X[[j]], type = "+", a = a, b = b)
  for (j in which(index == "-")) X[[j]] = rescale(X[[j]], type = "-", a = a, b = b)
  X
}
