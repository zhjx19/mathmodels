# ============================================================================
# Multiple-imputation inference: Rubin pooling and a convenience wrapper
# for lm/glm fits. Pool point estimates and standard errors across
# imputations according to Rubin's rules (with the Barnard-Rubin
# small-sample df).
# ============================================================================

#' Pool Estimates via Rubin's Rules
#'
#' Pools point estimates and standard errors from multiple imputations
#' according to Rubin's rules. Degrees of freedom follow the Barnard-Rubin
#' (1999) small-sample adjustment when complete-data degrees of freedom are
#' supplied via `edf`.
#'
#' @param est Numeric vector of point estimates, one per imputation.
#' @param err Numeric vector of corresponding standard errors.
#' @param alpha Significance level for the confidence interval. Default 0.05.
#' @param edf Complete-data degrees of freedom (positive scalar or `Inf`).
#'   When finite, the Barnard-Rubin adjusted df are used; default `Inf`
#'   applies the Rubin large-sample df.
#'
#' @return A list of class `rubin_pool`: `estimate`, `std_error`, `m`,
#'   `var_between`, `var_within`, `var_total`, `df`, `r`, `fmi`, `rel_eff`,
#'   `lcl`, `ucl`, `t`, `p_value` (two-sided test against 0), `alpha`,
#'   `edf`.
#'
#' @details
#' With Q_i and U_i the estimate and variance from imputation i, the pooled
#' estimate is the mean of Q_i; total variance T = Ubar + (1 + 1/m) B
#' where Ubar is the mean variance and B the between-imputation variance.
#' The fraction of missing information uses the standard Rubin definition
#' (r + 2/(df+3))/(r+1); the relative efficiency is (1 + lambda/m)^-1.
#' The fraction of missing information and relative efficiency use the
#' final (Barnard-Rubin adjusted) degrees of freedom when `edf` is finite.
#'
#' @examples
#' est = c(0.345, 0.303, 0.298, 0.378)
#' err = c(0.056, 0.101, 0.099, 0.034)
#' rubin_pool(est, err)
#' rubin_pool(est, err, edf = 120)
#' @export
rubin_pool = function(est, err, alpha = 0.05, edf = Inf) {
  if (!is.numeric(est) || !is.numeric(err)) {
    stop("est and err must be numeric vectors.", call. = FALSE)
  }
  if (length(est) != length(err)) {
    stop("est and err must have the same length.", call. = FALSE)
  }
  if (length(est) < 1L) {
    stop("est and err must have at least one value.", call. = FALSE)
  }
  if (anyNA(est) || anyNA(err) || any(!is.finite(est)) ||
      any(!is.finite(err))) {
    stop("est and err must contain only finite, non-missing values.",
         call. = FALSE)
  }
  if (any(err < 0)) {
    stop("err must be non-negative.", call. = FALSE)
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single value in (0, 1).", call. = FALSE)
  }
  if (length(edf) != 1L || is.na(edf) || edf <= 0) {
    stop("edf must be a single positive value or Inf.", call. = FALSE)
  }

  m = length(est)
  estimate = mean(est)
  var_within = mean(err^2)
  var_between = if (m > 1L) sum((est - estimate)^2) / (m - 1) else 0
  var_total = var_within + var_between + var_between / m

  if (var_between == 0 && var_within == 0) {
    r = 0
    df = Inf
  } else {
    r = if (var_within == 0) Inf else
      (1 + 1 / m) * var_between / var_within
    df = if (var_between == 0) Inf else
      (m - 1) * (1 + var_within / (var_between + var_between / m))^2
  }

  gamma = if (var_total == 0) 0 else (1 + 1 / m) * var_between / var_total
  if (is.finite(edf)) {
    v_obs = (1 - gamma) * edf * (edf + 1) / (edf + 3)
    df = if (is.infinite(df)) v_obs else 1 / (1 / df + 1 / v_obs)
  }

  fmi = if (is.infinite(r)) {
    1
  } else if (is.infinite(df)) {
    r / (r + 1)
  } else {
    (r + 2 / (df + 3)) / (r + 1)
  }

  std_error = sqrt(var_total)
  crit = stats::qt(1 - alpha / 2, df = df)
  rel_eff = 1 / (1 + fmi / m)
  t_val = if (std_error == 0) {
    if (estimate == 0) NA_real_ else sign(estimate) * Inf
  } else {
    estimate / std_error
  }
  p_value = if (is.na(t_val)) NA_real_ else
    2 * stats::pt(abs(t_val), df = df, lower.tail = FALSE)

  res = list(
    estimate = estimate,
    std_error = std_error,
    m = m,
    var_between = var_between,
    var_within = var_within,
    var_total = var_total,
    df = df,
    r = r,
    fmi = fmi,
    rel_eff = rel_eff,
    lcl = estimate - crit * std_error,
    ucl = estimate + crit * std_error,
    t = t_val,
    p_value = p_value,
    alpha = alpha,
    edf = edf
  )
  class(res) = c("rubin_pool", "list")
  res
}

#' @export
print.rubin_pool = function(x, ...) {
  cat("Rubin pooling (m =", x$m, ")\n")
  cat(sprintf("  estimate = %.6f, std_error = %.6f, df = %.1f\n",
              x$estimate, x$std_error, x$df))
  cat(sprintf(
    "  %.1f%% CI: [%.6f, %.6f], t = %.4f, p = %.4g, FMI = %.4f\n",
    100 * (1 - x$alpha), x$lcl, x$ucl, x$t, x$p_value, x$fmi))
  invisible(x)
}

# ----------------------------------------------------------------------------
# Convenience pooling for lm/glm fits
# ----------------------------------------------------------------------------

#' Pool lm/glm Fits Across Imputations
#'
#' Extracts coefficients and standard errors from a list of `lm`/`glm`
#' fits (one per imputation, usually built with [impute_multiple()]) and
#' pools each coefficient via [rubin_pool()], automatically passing the
#' complete-data residual degrees of freedom so the Barnard-Rubin
#' small-sample adjustment applies.
#'
#' @param fits A list of `lm`/`glm` fits, same model on the imputed data.
#' @param alpha Significance level. Default 0.05.
#'
#' @return A tibble with one row per coefficient: `term, estimate,
#'   std_error, lcl, ucl, df, t, p_value, fmi, r, rel_eff`.
#'
#' @examples
#' set.seed(1)
#' d = data.frame(x = rnorm(40), z = rnorm(40))
#' d$y = 0.5 * d$x + rnorm(40)
#' d$x[c(4, 9)] = NA
#' imps = impute_multiple(d, .cols = x, m = 5, maxit = 3, seed = 2)
#' fits = lapply(imps, \(dat) lm(y ~ x, dat))
#' pool_fit(fits)
#' @export
pool_fit = function(fits, alpha = 0.05) {
  if (!is.list(fits) || length(fits) < 1L ||
      !all(vapply(fits, \(f) inherits(f, "lm"), logical(1)))) {
    stop("fits must be a non-empty list of lm/glm fits.", call. = FALSE)
  }
  coefs = lapply(fits, stats::coef)
  term_names = names(coefs[[1]])
  same = all(vapply(coefs[-1], \(cn) identical(names(cn), term_names),
                    logical(1)))
  if (!same) {
    stop("coefficient names differ across fits.", call. = FALSE)
  }
  edf = stats::df.residual(fits[[1]])
  if (!is.finite(edf)) edf = Inf

  # name-based extraction survives aliased fits (NA coefficients)
  est = vapply(fits, \(f) {
    unname(stats::coef(f)[term_names])
  }, numeric(length(term_names)))
  err = vapply(fits, \(f) {
    se = tryCatch(sqrt(diag(stats::vcov(f))),
                  error = function(e) rep(NA_real_, length(term_names)))
    unname(se[term_names])
  }, numeric(length(term_names)))
  # vapply simplifies to a plain m-vector when p = 1; orient to m x p
  est = if (is.matrix(est)) t(est) else matrix(est, ncol = 1)
  err = if (is.matrix(err)) t(err) else matrix(err, ncol = 1)

  rows = vector("list", length(term_names))
  dropped = character(0)
  for (j in seq_along(term_names)) {
    ok = !is.na(est[, j]) & !is.na(err[, j])
    if (!all(ok)) {
      dropped = c(dropped, term_names[j])
      next
    }
    pooled = rubin_pool(est[ok, j], err[ok, j], alpha = alpha, edf = edf)
    rows[[j]] = tibble::tibble(
      term = term_names[j],
      estimate = pooled$estimate,
      std_error = pooled$std_error,
      lcl = pooled$lcl,
      ucl = pooled$ucl,
      df = pooled$df,
      t = pooled$t,
      p_value = pooled$p_value,
      fmi = pooled$fmi,
      r = pooled$r,
      rel_eff = pooled$rel_eff
    )
  }
  if (length(dropped) > 0) {
    warning("NA coefficient(s) dropped: ", paste(dropped, collapse = ", "),
            ".", call. = FALSE)
  }
  out = dplyr::bind_rows(rows)
  class(out) = c("pool_fit", class(out))
  out
}

#' @export
print.pool_fit = function(x, ...) {
  cat("Pooled estimates (Rubin's rules)\n")
  y = x[, c("term", "estimate", "std_error", "lcl", "ucl", "p_value",
            "fmi")]
  class(y) = setdiff(class(y), "pool_fit")
  print(y)
  invisible(x)
}
