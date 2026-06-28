# =============================================================================
# ts_models.R
# Time Series Modeling Functions
# Dependencies: forecast, tseries, rugarch
# =============================================================================

# ---------------------------------------------------------------------------
# Internal Helpers
# ---------------------------------------------------------------------------

#' Validate and coerce input to ts object
#' @keywords internal
.validate_ts_input = function(x, frequency = NULL) {
  if (inherits(x, "ts")) {
    return(x)
  }
  if (is.numeric(x)) {
    freq = if (!is.null(frequency)) frequency else 1L
    return(stats::ts(x, frequency = freq))
  }
  if (is.data.frame(x) || tibble::is_tibble(x)) {
    if (ncol(x) < 1) stop("Data frame must have at least one numeric column.")
    vals = x[[ncol(x)]]
    freq = if (!is.null(frequency)) frequency else 1L
    return(stats::ts(vals, frequency = freq))
  }
  stop("'x' must be a numeric vector, ts object, or data frame.")
}

#' Build a tidy coefficient tibble from a named numeric vector
#' @keywords internal
.tidy_coef = function(coefs, ses = NULL, model_name = "") {
  tbl = tibble::tibble(
    model    = model_name,
    term     = names(coefs),
    estimate = as.numeric(coefs)
  )
  if (!is.null(ses)) tbl$std_error = as.numeric(ses)
  tbl
}

#' Tidy residuals tibble
#' @keywords internal
.tidy_residuals = function(resid_vec, model_name = "") {
  tibble::tibble(
    model    = model_name,
    index    = seq_along(resid_vec),
    residual = as.numeric(resid_vec)
  )
}

# =============================================================================
# S1. ts_transform()  — Preprocessing: Box-Cox, Log, Differencing
# =============================================================================

#' Transform a Time Series for Stationarity
#'
#' Applies variance-stabilising transformation (Box-Cox / log) and/or
#' differencing, returning the transformed series together with all parameters
#' needed to back-transform forecasts.
#'
#' @param x A numeric vector or \code{ts} object.
#' @param method Character. Variance transformation:
#'   \code{"none"} (default), \code{"log"}, \code{"boxcox"} (auto-selects
#'   lambda via \code{forecast::BoxCox.lambda()}).
#' @param lambda Numeric. Fixed Box-Cox lambda. Ignored unless
#'   \code{method = "boxcox"}. If \code{NULL} (default), lambda is estimated.
#' @param diff Integer. Number of regular differences (default \code{0}).
#' @param seasonal_diff Integer. Number of seasonal differences (default
#'   \code{0}).
#' @param frequency Integer. Series frequency, required when \code{x} is a
#'   plain numeric vector and \code{seasonal_diff > 0}.
#' @param verbose Logical. Print transformation summary (default \code{TRUE}).
#'
#' @return A named list:
#'   \describe{
#'     \item{transformed}{Transformed \code{ts} object ready for modelling.}
#'     \item{params}{Named list with \code{method}, \code{lambda},
#'       \code{diff}, \code{seasonal_diff}, \code{frequency} — everything
#'       needed by \code{ts_back_transform()}.}
#'     \item{summary}{One-row tibble describing each step applied.}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' res = ts_transform(AirPassengers, method = "log", diff = 1,
#'                     seasonal_diff = 1)
#' res$summary
#' plot_ts(res$transformed, title = "Transformed AirPassengers")
#'
#' @seealso \code{\link{ts_back_transform}}
#' @export
ts_transform = function(x,
                         method        = c("none", "log", "boxcox"),
                         lambda        = NULL,
                         diff          = 0L,
                         seasonal_diff = 0L,
                         frequency     = NULL,
                         verbose       = TRUE) {
  method = match.arg(method)
  x_ts   = .validate_ts_input(x, frequency)
  freq   = stats::frequency(x_ts)

  lam_used  = NA_real_
  x_out     = x_ts
  steps     = character(0)

  # --- Step 1: Variance stabilisation ---
  if (method == "log") {
    if (any(x_out <= 0, na.rm = TRUE))
      stop("Log transform requires all values > 0.")
    x_out    = log(x_out)
    steps    = c(steps, "log(x)")
  } else if (method == "boxcox") {
    lam_used = if (!is.null(lambda)) lambda else
      forecast::BoxCox.lambda(x_out)
    x_out    = forecast::BoxCox(x_out, lam_used)
    steps    = c(steps, sprintf("BoxCox(lambda=%.4f)", lam_used))
  }

  # Save variance-stabilised-only series (before differencing)
  # for ts_back_transform() to use in un-differencing
  x_vs = x_out

  # --- Step 2: Seasonal differencing ---
  if (seasonal_diff > 0L) {
    for (i in seq_len(seasonal_diff)) {
      x_out = diff(x_out, lag = freq)
      steps = c(steps, sprintf("seasonal_diff(lag=%d)", freq))
    }
  }

  # --- Step 3: Regular differencing ---
  if (diff > 0L) {
    for (i in seq_len(diff)) {
      x_out = diff(x_out)
      steps = c(steps, "diff(lag=1)")
    }
  }
  params = list(
    method        = method,
    lambda        = lam_used,
    diff          = diff,
    seasonal_diff = seasonal_diff,
    frequency     = freq,
    n_original    = length(x_ts),
    x_vs          = x_vs
  )

  summary_tbl = tibble::tibble(
    original_length    = length(x_ts),
    transformed_length = length(x_out),
    steps_applied      = if (length(steps) == 0) "none" else paste(steps, collapse = " -> "),
    lambda             = lam_used
  )

  if (verbose) {
    cat("\n=== ts_transform() ===\n")
    cat("Steps:", summary_tbl$steps_applied, "\n")
    cat("Length:", summary_tbl$original_length, "->",
        summary_tbl$transformed_length, "\n")
    if (method == "boxcox") cat("Lambda:", round(lam_used, 4), "\n")
  }

  list(transformed = x_out, params = params, summary = summary_tbl)
}


#' Back-Transform Forecasts to the Original Scale
#'
#' Reverses the operations applied by \code{ts_transform()} to restore
#' point forecasts and confidence intervals to the original scale.
#' Differencing is un-done by cumulative summation using the tail of the
#' original series as the starting level.
#'
#' @param forecast_tbl Tibble from \code{ts_forecast()} (columns: \code{step},
#'   \code{forecast}, optionally \code{lo_*} and \code{hi_*}).
#' @param params The \code{$params} element from a \code{ts_transform()} result.
#' @param x_original The original (pre-transform) \code{ts} or numeric vector,
#'   needed to supply the last observed level(s) for un-differencing.
#'
#' @return The \code{forecast_tbl} with \code{forecast} and all interval
#'   columns converted to the original scale.
#'
#' @examples
#' data(AirPassengers)
#' tr  = ts_transform(AirPassengers, method = "log", diff = 1,
#'                     seasonal_diff = 1, verbose = FALSE)
#' fit = ts_sarima(tr$transformed)
#' fc  = ts_forecast(fit, h = 12)
#' fc_orig = ts_back_transform(fc, tr$params, AirPassengers)
#'
#' @export
ts_back_transform = function(forecast_tbl, params, x_original) {
  x_orig = as.numeric(x_original)
  freq   = params$frequency

  # Columns to invert (all numeric except 'step')
  num_cols = setdiff(names(forecast_tbl)[sapply(forecast_tbl, is.numeric)],
                      "step")

  tbl = forecast_tbl

  # --- Un-difference (regular) ---
  # Uses variance-stabilised level (params$x_vs), not raw scale
  if (params$diff > 0L) {
    last_val = utils::tail(params$x_vs, 1)
    for (col in num_cols) {
      tbl[[col]] = cumsum(c(last_val, tbl[[col]]))[-1]
    }
  }

  # --- Un-seasonal-difference ---
  if (params$seasonal_diff > 0L) {
    last_season = utils::tail(params$x_vs, freq)
    for (col in num_cols) {
      out = numeric(nrow(tbl))
      for (i in seq_len(nrow(tbl))) {
        base_val = if (i <= freq) last_season[i] else out[i - freq]
        out[i]   = tbl[[col]][i] + base_val
      }
      tbl[[col]] = out
    }
  }

  # --- Un-transform variance stabilisation ---
  if (params$method == "log") {
    for (col in num_cols) tbl[[col]] = exp(tbl[[col]])
  } else if (params$method == "boxcox") {
    for (col in num_cols)
      tbl[[col]] = forecast::InvBoxCox(tbl[[col]], params$lambda)
  }

  tbl
}


# =============================================================================
# 1. ts_test()  — Stationarity Tests
# =============================================================================

#' Stationarity Tests for a Time Series
#'
#' Runs ADF, KPSS, and PP tests and returns a tidy summary data frame.
#'
#' @param x A numeric vector, \code{ts} object, or single-column data frame.
#' @param adf_lags Integer. Max lag for ADF (default: \code{NULL} = auto via
#'   \code{trunc((length(x)-1)^(1/3))}).
#' @param kpss_type Character. \code{"Level"} (default) or \code{"Trend"}
#'   for KPSS null hypothesis.
#' @param pp_type  Character. \code{"Z(t_alpha)"} (default) or \code{"Z(alpha)"}.
#' @param verbose Logical. Print test results to console (default \code{TRUE}).
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{test}{Test name (ADF / KPSS / PP)}
#'     \item{null_hypothesis}{Plain-English description of H0}
#'     \item{statistic}{Test statistic}
#'     \item{p_value}{p-value (or NA when only bounds are available)}
#'     \item{conclusion}{Character: "Stationary" or "Non-stationary"}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' ts_test(AirPassengers)
#'
#' @export
ts_test = function(x,
                   adf_lags  = NULL,
                   kpss_type = c("Level", "Trend"),
                   pp_type   = c("Z(t_alpha)", "Z(alpha)"),
                   verbose   = TRUE) {
  pp_type   = match.arg(pp_type)

  x_ts = .validate_ts_input(x)
  n    = length(x_ts)

  # --- ADF ---
  adf_lag = if (is.null(adf_lags)) trunc((n - 1)^(1 / 3)) else adf_lags
  adf_res = tseries::adf.test(x_ts, k = adf_lag)
  adf_row = tibble::tibble(
    test           = "ADF",
    null_hypothesis = "Unit root present (non-stationary)",
    statistic      = unname(adf_res$statistic),
    p_value        = adf_res$p.value,
    conclusion     = ifelse(adf_res$p.value < 0.05,
                            "Stationary", "Non-stationary")
  )

  # --- KPSS ---
  kpss_res = tseries::kpss.test(x_ts, null = kpss_type)
  kpss_row = tibble::tibble(
    test           = "KPSS",
    null_hypothesis = "Series is stationary",
    statistic      = unname(kpss_res$statistic),
    p_value        = kpss_res$p.value,
    conclusion     = ifelse(kpss_res$p.value >= 0.05,
                            "Stationary", "Non-stationary")
  )

  # --- PP ---
  pp_res = tseries::pp.test(x_ts, type = pp_type)
  pp_row = tibble::tibble(
    test           = "PP",
    null_hypothesis = "Unit root present (non-stationary)",
    statistic      = unname(pp_res$statistic),
    p_value        = pp_res$p.value,
    conclusion     = ifelse(pp_res$p.value < 0.05,
                            "Stationary", "Non-stationary")
  )

  result = dplyr::bind_rows(adf_row, kpss_row, pp_row)

  if (verbose) {
    cat("\n=== Stationarity Tests ===\n")
    print(result, n = Inf)
    cat("\nNote: All tests use alpha = 0.05.\n")
  }

  invisible(result)
}


# =============================================================================
# 2. ts_stl()  — STL Decomposition
# =============================================================================

#' STL (Seasonal + Trend + Loess) Decomposition
#'
#' Wraps \code{stats::stl()} and returns both a tidy long-format tibble and
#' the raw \code{stl} object for further use.
#'
#' @param x A \code{ts} object with \code{frequency > 1}, or a numeric vector
#'   with \code{frequency} specified.
#' @param frequency Integer. Required if \code{x} is a plain numeric vector.
#' @param s_window Seasonal smoothing window passed to \code{stl()}.
#'   Use \code{"periodic"} (default) for purely periodic seasonality.
#' @param robust Logical. Use robust fitting in STL (default \code{TRUE}).
#' @param ... Additional arguments forwarded to \code{stats::stl()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{components}{Tidy \code{tibble} with columns
#'       \code{index, observed, trend, seasonal, remainder}.}
#'     \item{strength}{Tibble: \code{seasonal_strength} and
#'       \code{trend_strength} (Wang et al. 2006 metrics).}
#'     \item{model}{Raw \code{stl} object.}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' res = ts_stl(AirPassengers)
#' res$components
#'
#' @export
ts_stl = function(x,
                   frequency = NULL,
                   s_window  = "periodic",
                   robust    = TRUE,
                   ...) {
  x_ts = .validate_ts_input(x, frequency)
  if (stats::frequency(x_ts) <= 1)
    stop("STL decomposition requires frequency > 1 (seasonal data).")

  stl_fit = stats::stl(x_ts, s.window = s_window, robust = robust, ...)
  comp    = stl_fit$time.series  # matrix: seasonal, trend, remainder

  components = tibble::tibble(
    index     = seq_len(nrow(comp)),
    observed  = as.numeric(x_ts),
    trend     = as.numeric(comp[, "trend"]),
    seasonal  = as.numeric(comp[, "seasonal"]),
    remainder = as.numeric(comp[, "remainder"])
  )

  # Strength metrics (Hyndman & Athanasopoulos §2.3 / Wang et al. 2006)
  var_rem  = stats::var(components$remainder, na.rm = TRUE)
  s_strength = max(0, 1 - var_rem /
                      stats::var(components$seasonal + components$remainder,
                                 na.rm = TRUE))
  t_strength = max(0, 1 - var_rem /
                      stats::var(components$trend + components$remainder,
                                 na.rm = TRUE))

  strength = tibble::tibble(
    seasonal_strength = round(s_strength, 4),
    trend_strength    = round(t_strength, 4)
  )

  list(components = components, strength = strength, model = stl_fit)
}


# =============================================================================
# 3. ts_ets()  — Exponential Smoothing (ETS)
# =============================================================================

#' ETS (Error, Trend, Seasonality) Exponential Smoothing
#'
#' Wraps \code{forecast::ets()} and returns a tidy result list.
#'
#' @param x A numeric vector or \code{ts} object.
#' @param model Character ETS model string, e.g. \code{"ZZZ"} (auto-select,
#'   default), \code{"AAN"}, \code{"AAA"}.
#' @param frequency Integer. Required when \code{x} is a plain numeric vector.
#' @param ... Additional arguments forwarded to \code{forecast::ets()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{model_info}{One-row tibble: model type, AIC, AICc, BIC,
#'       log-likelihood.}
#'     \item{parameters}{Tidy tibble of smoothing parameters and initial
#'       states.}
#'     \item{fitted}{Tibble with \code{index, observed, fitted, residual}.}
#'     \item{model}{Raw \code{ets} object.}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' res = ts_ets(AirPassengers)
#' res$model_info
#'
#' @export
ts_ets = function(x, model = "ZZZ", frequency = NULL, ...) {
  x_ts    = .validate_ts_input(x, frequency)
  ets_fit = forecast::ets(x_ts, model = model, ...)

  model_info = tibble::tibble(
    model_type  = ets_fit$method,
    log_lik     = round(ets_fit$loglik, 4),
    aic         = round(ets_fit$aic, 4),
    aicc        = round(ets_fit$aicc, 4),
    bic         = round(ets_fit$bic, 4)
  )

  # Smoothing parameters
  pars = c(ets_fit$par)
  parameters = tibble::tibble(
    term     = names(pars),
    estimate = round(as.numeric(pars), 6)
  )

  fitted_tbl = tibble::tibble(
    index    = seq_along(x_ts),
    observed = as.numeric(x_ts),
    fitted   = as.numeric(ets_fit$fitted),
    residual = as.numeric(ets_fit$residuals)
  )

  list(model_info = model_info,
       parameters = parameters,
       fitted     = fitted_tbl,
       model      = ets_fit)
}


# =============================================================================
# 4. ts_sarima()  — SARIMA Mean Modeling
# =============================================================================

#' SARIMA Model Fitting
#'
#' Auto-selects or fits a user-specified SARIMA model via \code{forecast::auto.arima()}
#' or \code{forecast::Arima()}, returning a comprehensive tidy result.
#'
#' @param x A numeric vector or \code{ts} object.
#' @param order Integer vector of length 3: \code{c(p, d, q)}.
#'   \code{NULL} (default) triggers automatic selection.
#' @param seasonal A list \code{list(order = c(P,D,Q), period = m)} or
#'   \code{NULL} for automatic selection.
#' @param frequency Integer. Required when \code{x} is a plain vector.
#' @param stepwise Logical. Use stepwise search in auto.arima (default \code{TRUE}).
#' @param approximation Logical. Use approximation in auto.arima (default \code{TRUE}).
#' @param ... Additional arguments forwarded to \code{forecast::auto.arima()}
#'   or \code{forecast::Arima()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{model_info}{One-row tibble: ARIMA order string, AIC, AICc, BIC,
#'       log-likelihood, sigma².}
#'     \item{coefficients}{Tidy tibble with estimate and std error.}
#'     \item{fitted}{Tibble: \code{index, observed, fitted, residual}.}
#'     \item{diagnostics}{Ljung-Box test result tibble.}
#'     \item{model}{Raw \code{Arima} object.}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' res = ts_sarima(log(AirPassengers))
#' res$model_info
#'
#' @export
ts_sarima = function(x,
                      order        = NULL,
                      seasonal     = NULL,
                      frequency    = NULL,
                      stepwise     = TRUE,
                      approximation = TRUE,
                      ...) {
  x_ts = .validate_ts_input(x, frequency)

  if (is.null(order)) {
    fit = forecast::auto.arima(x_ts,
                                stepwise     = stepwise,
                                approximation = approximation,
                                ...)
  } else {
    seas = if (is.null(seasonal)) list(order = c(0, 0, 0)) else seasonal
    fit  = forecast::Arima(x_ts, order = order, seasonal = seas, ...)
  }

  # Model info
  arima_str = paste0("ARIMA", paste0(fit$arma[c(1,6,2)], collapse=","),
                      if (any(fit$arma[c(3,7,4)] > 0))
                        paste0("(", paste0(fit$arma[c(3,7,4)], collapse=","),
                               ")[", fit$arma[5], "]") else "")

  model_info = tibble::tibble(
    model_type = arima_str,
    log_lik    = round(fit$loglik, 4),
    aic        = round(fit$aic, 4),
    aicc       = round(fit$aicc, 4),
    bic        = round(fit$bic, 4),
    sigma2     = round(fit$sigma2, 6)
  )

  # Coefficients
  coef_vec = stats::coef(fit)
  se_vec   = sqrt(diag(fit$var.coef))
  coefficients = tibble::tibble(
    term      = names(coef_vec),
    estimate  = round(as.numeric(coef_vec), 6),
    std_error = round(as.numeric(se_vec), 6),
    t_stat    = round(as.numeric(coef_vec) / as.numeric(se_vec), 4)
  )

  # Fitted values
  fitted_tbl = tibble::tibble(
    index    = seq_along(x_ts),
    observed = as.numeric(x_ts),
    fitted   = as.numeric(fit$fitted),
    residual = as.numeric(fit$residuals)
  )

  # Ljung-Box diagnostics
  lb = stats::Box.test(fit$residuals, lag = min(20, length(x_ts) / 5),
                        type = "Ljung-Box", fitdf = length(coef_vec))
  diagnostics = tibble::tibble(
    test           = "Ljung-Box",
    lag            = lb$parameter,
    statistic      = round(lb$statistic, 4),
    p_value        = round(lb$p.value, 4),
    conclusion     = ifelse(lb$p.value > 0.05,
                            "No autocorrelation (residuals OK)",
                            "Autocorrelation detected")
  )

  list(model_info    = model_info,
       coefficients  = coefficients,
       fitted        = fitted_tbl,
       diagnostics   = diagnostics,
       model         = fit)
}


# =============================================================================
# 5. ts_garch()  — Pure GARCH Variance Modeling
# =============================================================================

#' GARCH Variance Modeling
#'
#' Fits a GARCH(p,q) model (optionally with ARMA mean) using
#' \code{rugarch::ugarchfit()} and returns tidy results.
#'
#' @param x A numeric vector or \code{ts} object (typically returns/residuals).
#' @param garch_order Integer vector \code{c(p, q)}: GARCH order (default
#'   \code{c(1, 1)}).
#' @param arma_order Integer vector \code{c(p, q)}: ARMA order for the mean
#'   equation (default \code{c(0, 0)} = constant mean).
#' @param dist Character. Innovation distribution: \code{"norm"} (default),
#'   \code{"std"} (Student-t), \code{"ged"}, \code{"snorm"}, \code{"sstd"}.
#' @param frequency Integer. Required when \code{x} is a plain vector.
#' @param ... Additional arguments forwarded to \code{rugarch::ugarchspec()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{model_info}{One-row tibble: spec, distribution, log-likelihood,
#'       AIC, BIC.}
#'     \item{coefficients}{Tidy coefficient tibble with estimates, std errors,
#'       t-stats, p-values.}
#'     \item{fitted}{Tibble: \code{index, observed, sigma, variance,
#'       std_residual}.}
#'     \item{diagnostics}{ARCH-LM and Ljung-Box tests on squared residuals.}
#'     \item{model}{Raw \code{uGARCHfit} object.}
#'   }
#'
#' @examples
#' set.seed(42)
#' r = diff(log(AirPassengers))
#' res = ts_garch(r)
#' res$model_info
#'
#' @export
ts_garch = function(x,
                     garch_order = c(1, 1),
                     arma_order  = c(0, 0),
                     dist        = "norm",
                     frequency   = NULL,
                     ...) {
  if (!requireNamespace("rugarch", quietly = TRUE))
    stop("Package 'rugarch' is required for ts_garch(). Install it with install.packages('rugarch').")

  x_vec = as.numeric(.validate_ts_input(x, frequency))

  spec = rugarch::ugarchspec(
    variance.model = list(model = "sGARCH",
                          garchOrder = garch_order),
    mean.model     = list(armaOrder = arma_order,
                          include.mean = TRUE),
    distribution.model = dist,
    ...
  )

  fit = rugarch::ugarchfit(spec = spec, data = x_vec, solver = "hybrid")

  # Model info
  ic = rugarch::infocriteria(fit)
  model_info = tibble::tibble(
    model_type   = paste0("GARCH(", garch_order[1], ",", garch_order[2], ")"),
    distribution = dist,
    log_lik      = round(rugarch::likelihood(fit), 4),
    aic          = round(ic["Akaike", ], 4),
    bic          = round(ic["Bayes", ], 4)
  )

  # Coefficients
  coef_mat = rugarch::coef(fit)
  se_mat   = fit@fit$matcoef[, 2]
  tstat    = fit@fit$matcoef[, 3]
  pval     = fit@fit$matcoef[, 4]

  coefficients = tibble::tibble(
    term      = names(coef_mat),
    estimate  = round(as.numeric(coef_mat), 6),
    std_error = round(as.numeric(se_mat), 6),
    t_stat    = round(as.numeric(tstat), 4),
    p_value   = round(as.numeric(pval), 4)
  )

  # Fitted: sigma and std residuals
  sigma_vec  = as.numeric(rugarch::sigma(fit))
  z_vec      = as.numeric(rugarch::residuals(fit, standardize = TRUE))
  res_vec    = as.numeric(rugarch::residuals(fit))

  fitted_tbl = tibble::tibble(
    index        = seq_along(x_vec),
    observed     = x_vec,
    sigma        = sigma_vec,
    variance     = sigma_vec^2,
    residual     = res_vec,
    std_residual = z_vec
  )

  # Diagnostics: manual ARCH-LM (rugarch >= 1.5 removed ArchTest)
  n_obs   = length(z_vec)
  z2      = z_vec^2
  z2_lag  = stats::embed(z2, 11)[, -1]        # 10-lag matrix
  z2_dep  = z2[11:n_obs]

  lm_fit  = stats::lm(z2_dep ~ z2_lag)
  r2      = summary(lm_fit)$r.squared
  arch_stat = (n_obs - 10) * r2
  arch_pv   = stats::pchisq(arch_stat, df = 10, lower.tail = FALSE)

  lb_test = stats::Box.test(z2, lag = 10, type = "Ljung-Box")

  diagnostics = tibble::tibble(
    test      = c("ARCH-LM (sq.std.resid)", "Ljung-Box (sq.std.resid)"),
    lag       = c(10L, 10L),
    statistic = round(c(arch_stat, lb_test$statistic), 4),
    p_value   = round(c(arch_pv, lb_test$p.value), 4),
    conclusion = ifelse(c(arch_pv, lb_test$p.value) > 0.05,
                        "No ARCH effect remaining",
                        "ARCH effect still present")
  )

  list(model_info   = model_info,
       coefficients = coefficients,
       fitted       = fitted_tbl,
       diagnostics  = diagnostics,
       model        = fit)
}


# =============================================================================
# 6. ts_sarima_garch()  — Two-Stage SARIMA + GARCH
# =============================================================================

#' Two-Stage SARIMA-GARCH Joint Model
#'
#' Stage 1: fits a SARIMA model on the level series via \code{ts_sarima()}.
#' Stage 2: fits a GARCH model on the SARIMA residuals via \code{ts_garch()}.
#' Returns both sub-models plus a combined fitted tibble.
#'
#' @param x A numeric vector or \code{ts} object.
#' @param sarima_order Integer vector \code{c(p,d,q)} or \code{NULL} (auto).
#' @param sarima_seasonal List \code{list(order=c(P,D,Q), period=m)} or \code{NULL}.
#' @param garch_order Integer vector \code{c(p,q)} (default \code{c(1,1)}).
#' @param garch_dist Character. Innovation distribution for GARCH (default
#'   \code{"std"}).
#' @param frequency Integer. Required for plain numeric \code{x}.
#' @param ... Additional arguments forwarded to \code{ts_sarima()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{mean_model}{Result list from \code{ts_sarima()}.}
#'     \item{variance_model}{Result list from \code{ts_garch()}.}
#'     \item{fitted}{Combined tibble: \code{index, observed, mean_fitted,
#'       sigma, variance, std_residual}.}
#'     \item{model_info}{Combined one-row summary tibble.}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' res = ts_sarima_garch(log(AirPassengers))
#' res$model_info
#'
#' @export
ts_sarima_garch = function(x,
                            sarima_order    = NULL,
                            sarima_seasonal = NULL,
                            garch_order     = c(1, 1),
                            garch_dist      = "std",
                            frequency       = NULL,
                            ...) {
  x_ts = .validate_ts_input(x, frequency)

  # Stage 1: SARIMA on level series
  mean_res = ts_sarima(x_ts,
                        order    = sarima_order,
                        seasonal = sarima_seasonal,
                        ...)

  sarima_resid = mean_res$fitted$residual

  # Stage 2: GARCH on SARIMA residuals
  var_res = ts_garch(sarima_resid,
                      garch_order = garch_order,
                      dist        = garch_dist)

  # Combined fitted tibble
  combined = tibble::tibble(
    index        = seq_along(x_ts),
    observed     = as.numeric(x_ts),
    mean_fitted  = mean_res$fitted$fitted,
    sigma        = var_res$fitted$sigma,
    variance     = var_res$fitted$variance,
    std_residual = var_res$fitted$std_residual
  )

  model_info = tibble::tibble(
    mean_model     = mean_res$model_info$model_type,
    variance_model = var_res$model_info$model_type,
    garch_dist     = garch_dist,
    sarima_aic     = mean_res$model_info$aic,
    garch_aic      = var_res$model_info$aic
  )

  list(mean_model     = mean_res,
       variance_model = var_res,
       fitted         = combined,
       model_info     = model_info)
}


# =============================================================================
# 7. ts_forecast()  — Generic Forecast Function
# =============================================================================

#' Generate Forecasts from a Fitted Time Series Model
#'
#' A unified interface for forecasting from ETS, SARIMA, GARCH, or
#' SARIMA-GARCH model results produced by the \code{ts_*()} functions.
#'
#' @param model_result A result list from \code{ts_ets()}, \code{ts_sarima()},
#'   \code{ts_garch()}, or \code{ts_sarima_garch()}.
#' @param h Integer. Forecast horizon (number of steps ahead).
#' @param level Numeric vector. Confidence levels in percent (default
#'   \code{c(80, 95)}).
#' @param x_ts Original \code{ts} object (required for GARCH/SARIMA-GARCH to
#'   preserve time attributes). If \code{NULL}, integer steps are used.
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{step}{Forecast step (1 to h).}
#'     \item{forecast}{Point forecast.}
#'     \item{lo_<level>}{Lower bound for each confidence level.}
#'     \item{hi_<level>}{Upper bound for each confidence level.}
#'     \item{sigma}{Predicted conditional std dev (GARCH models only).}
#'   }
#'
#' @examples
#' data(AirPassengers)
#' fit = ts_sarima(log(AirPassengers))
#' ts_forecast(fit, h = 24)
#'
#' @export
ts_forecast = function(model_result, h = 12, level = c(80, 95), x_ts = NULL) {

  # ---- Detect model type ----
  has_slot = function(obj, nm) !is.null(obj[[nm]])

  # SARIMA-GARCH (has both mean_model and variance_model)
  if (has_slot(model_result, "mean_model") &&
      has_slot(model_result, "variance_model")) {
    return(.forecast_sarima_garch(model_result, h, level))
  }

  # GARCH only (has $model of class uGARCHfit)
  if (inherits(model_result$model, "uGARCHfit")) {
    return(.forecast_garch(model_result, h, level))
  }

  # ETS or SARIMA (forecast objects)
  if (inherits(model_result$model, c("ets", "ARIMA", "Arima"))) {
    return(.forecast_ets_sarima(model_result, h, level))
  }

  stop("Unrecognised model_result type. Use output from ts_ets(), ts_sarima(), ts_garch(), or ts_sarima_garch().")
}

#' @keywords internal
.forecast_ets_sarima = function(res, h, level) {
  fc   = forecast::forecast(res$model, h = h, level = level)
  tbl  = tibble::tibble(step = seq_len(h),
                         forecast = as.numeric(fc$mean))
  for (lv in level) {
    idx = which(fc$level == lv)
    tbl[[paste0("lo_", lv)]] = as.numeric(fc$lower[, idx])
    tbl[[paste0("hi_", lv)]] = as.numeric(fc$upper[, idx])
  }
  tbl
}

#' @keywords internal
.forecast_garch = function(res, h, level) {
  if (!requireNamespace("rugarch", quietly = TRUE))
    stop("Package 'rugarch' required.")
  fc   = rugarch::ugarchforecast(res$model, n.ahead = h)
  mu   = as.numeric(rugarch::fitted(fc))
  sig  = as.numeric(rugarch::sigma(fc))

  tbl = tibble::tibble(step = seq_len(h), forecast = mu, sigma = sig)
  for (lv in level) {
    z = stats::qnorm(0.5 + lv / 200)
    tbl[[paste0("lo_", lv)]] = mu - z * sig
    tbl[[paste0("hi_", lv)]] = mu + z * sig
  }
  tbl
}

#' @keywords internal
.forecast_sarima_garch = function(res, h, level) {
  mean_fc = .forecast_ets_sarima(res$mean_model, h, level)
  var_fc  = .forecast_garch(res$variance_model, h, level)
  # Replace intervals with GARCH-based intervals around SARIMA mean
  tbl = tibble::tibble(step = seq_len(h),
                        forecast = mean_fc$forecast,
                        sigma    = var_fc$sigma)
  for (lv in level) {
    z = stats::qnorm(0.5 + lv / 200)
    tbl[[paste0("lo_", lv)]] = mean_fc$forecast - z * var_fc$sigma
    tbl[[paste0("hi_", lv)]] = mean_fc$forecast + z * var_fc$sigma
  }
  tbl
}

# =============================================================================
# ts_plots.R
# Time Series Visualization Functions for mathmodels Package
# Dependencies: ggplot2, patchwork, dplyr, tibble
# =============================================================================

# Silence R CMD check NOTE for computed ggplot2 aesthetics
utils::globalVariables("density")

# ---------------------------------------------------------------------------
# Internal Helpers
# ---------------------------------------------------------------------------

#' Shared ggplot2 theme for mathmodels time series plots
#' @keywords internal
.ts_theme = function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = base_size - 1),
      axis.title       = ggplot2::element_text(size = base_size - 1),
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}

#' Convert ts/numeric to tidy tibble with index
#' @keywords internal
.ts_to_tbl = function(x) {
  x_vec = as.numeric(x)
  tibble::tibble(index = seq_along(x_vec), value = x_vec)
}

# Shared colour palette
.PALETTE = list(
  observed  = "#2C3E50",
  fitted    = "#2980B9",
  forecast  = "#E74C3C",
  interval  = "#E74C3C",
  sigma     = "#8E44AD",
  trend     = "#27AE60",
  seasonal  = "#F39C12",
  remainder = "#95A5A6",
  ref_line  = "#BDC3C7"
)


# =============================================================================
# 1. plot_ts()  — Basic Time Series Line Plot
# =============================================================================

#' Plot a Time Series
#'
#' Renders a clean line chart for one or more time series.
#'
#' @param x A numeric vector, \code{ts} object, tidy tibble with columns
#'   \code{(index, value)}, or a named list of such objects for multi-series.
#' @param title Character. Plot title.
#' @param subtitle Character. Plot subtitle.
#' @param x_lab Character. x-axis label (default \code{"Time"}).
#' @param y_lab Character. y-axis label (default \code{"Value"}).
#' @param colour Character. Line colour (ignored for multi-series).
#' @param add_points Logical. Add point markers (default \code{FALSE}).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(AirPassengers)
#' plot_ts(AirPassengers, title = "Air Passengers", y_lab = "Thousands")
#'
#' @export
plot_ts = function(x,
                    title      = "Time Series",
                    subtitle   = NULL,
                    x_lab      = "Time",
                    y_lab      = "Value",
                    colour     = .PALETTE$observed,
                    add_points = FALSE) {

  # Multi-series (named list) branch
  if (is.list(x) && !is.data.frame(x)) {
    tbl = purrr::imap_dfr(x, function(s, nm) {
      d = .ts_to_tbl(s); d$series = nm; d
    })
    p = ggplot2::ggplot(tbl, ggplot2::aes(x = .data[["index"]], y = .data[["value"]],
                                           colour = .data[["series"]])) +
      ggplot2::geom_line(linewidth = 0.8)
    if (add_points) p = p + ggplot2::geom_point(size = 1.5)
    p = p + ggplot2::scale_colour_brewer(palette = "Set2",
                                          name = "Series")
  } else {
    tbl = if (is.data.frame(x) && "index" %in% names(x)) x else .ts_to_tbl(x)
    p = ggplot2::ggplot(tbl, ggplot2::aes(x = .data[["index"]], y = .data[["value"]])) +
      ggplot2::geom_line(colour = colour, linewidth = 0.9)
    if (add_points)
      p = p + ggplot2::geom_point(colour = colour, size = 1.5)
  }

  p + ggplot2::labs(title = title, subtitle = subtitle,
                    x = x_lab, y = y_lab) +
    .ts_theme()
}


# =============================================================================
# 2. plot_ts_forecast()  — General Forecast Plot
# =============================================================================

#' Plot Historical Series + Forecast with Confidence Bands
#'
#' Works with the tibble returned by \code{ts_forecast()}.
#'
#' @param observed A numeric vector or \code{ts} object of the historical data.
#' @param forecast_tbl Tibble from \code{ts_forecast()}.
#' @param level Numeric vector. Which confidence levels to shade (must match
#'   columns in \code{forecast_tbl}).
#' @param title Character. Plot title.
#' @param y_lab Character. y-axis label.
#' @param trim_n Integer. Show only the last \code{trim_n} observations from
#'   history (default \code{NULL} = show all).
#' @importFrom rlang .data
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(AirPassengers)
#' fit = ts_sarima(log(AirPassengers))
#' fc  = ts_forecast(fit, h = 24)
#' plot_ts_forecast(log(AirPassengers), fc, title = "Log Air Passengers Forecast")
#'
#' @export
plot_ts_forecast = function(observed,
                             forecast_tbl,
                             level    = c(80, 95),
                             title    = "Forecast",
                             y_lab    = "Value",
                             trim_n   = NULL) {

  obs_vec = as.numeric(observed)
  n_obs   = length(obs_vec)

  if (!is.null(trim_n)) {
    obs_vec = utils::tail(obs_vec, trim_n)
    n_obs   = length(obs_vec)
  }

  obs_tbl = tibble::tibble(index = seq_len(n_obs), value = obs_vec)
  n_start = n_obs + 1L

  fc_tbl = dplyr::mutate(forecast_tbl,
                          index = .data$step + n_obs - 1L + 1L)

  # Shade confidence bands (largest to smallest)
  level_sorted = sort(level, decreasing = TRUE)
  alpha_vals   = seq(0.12, 0.25, length.out = length(level_sorted))

  p = ggplot2::ggplot() +
    ggplot2::geom_line(data = obs_tbl,
                       ggplot2::aes(x = .data$index, y = .data$value),
                       colour = .PALETTE$observed, linewidth = 0.9)

  for (i in seq_along(level_sorted)) {
    lv  = level_sorted[i]
    lo_col = paste0("lo_", lv)
    hi_col = paste0("hi_", lv)
    if (lo_col %in% names(fc_tbl)) {
      p = p + ggplot2::geom_ribbon(
        data    = fc_tbl,
        ggplot2::aes(x = .data$index,
                     ymin = .data[[lo_col]],
                     ymax = .data[[hi_col]]),
        fill  = .PALETTE$interval,
        alpha = alpha_vals[i]
      )
    }
  }

  p + ggplot2::geom_line(data = fc_tbl,
                         ggplot2::aes(x = .data$index, y = .data$forecast),
                         colour = .PALETTE$forecast, linewidth = 1,
                         linetype = "dashed") +
    ggplot2::geom_vline(xintercept = n_obs + 0.5,
                        colour = .PALETTE$ref_line, linetype = "dotted") +
    ggplot2::labs(title    = title,
                  subtitle = paste0(nrow(forecast_tbl), "-step forecast | ",
                                    paste(level, collapse = "%, "), "% CI"),
                  x = "Time", y = y_lab) +
    .ts_theme()
}


# =============================================================================
# 3. plot_ts_sarima_garch()  — Dual-Panel Mean + Volatility
# =============================================================================

#' Dual-Axis Plot for SARIMA-GARCH: Mean + Conditional Volatility
#'
#' Upper panel shows observed values with mean fitted line; lower panel shows
#' the conditional standard deviation from the GARCH component.
#'
#' @param sg_result Result list from \code{ts_sarima_garch()}.
#' @param title Character. Overall plot title.
#' @param y_lab Character. y-axis label for upper panel.
#' @importFrom rlang .data
#'
#' @return A \code{patchwork} composite \code{ggplot} object.
#'
#' @examples
#' data(AirPassengers)
#' res = ts_sarima_garch(log(AirPassengers))
#' plot_ts_sarima_garch(res)
#'
#' @export
plot_ts_sarima_garch = function(sg_result,
                                 title = "SARIMA-GARCH: Mean & Volatility",
                                 y_lab = "Value") {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required. Install with install.packages('patchwork').")

  tbl = sg_result$fitted

  p_mean = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$observed),
                       colour = .PALETTE$observed, linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean_fitted),
                       colour = .PALETTE$fitted, linewidth = 0.9,
                       linetype = "dashed") +
    ggplot2::labs(subtitle = "Observed vs SARIMA Fitted",
                  x = NULL, y = y_lab) +
    .ts_theme()

  p_vol = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index, y = .data$sigma)) +
    ggplot2::geom_line(colour = .PALETTE$sigma, linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = mean(tbl$sigma),
                        colour = .PALETTE$ref_line, linetype = "dashed") +
    ggplot2::labs(subtitle = "Conditional Std Dev (GARCH sigma)",
                  x = "Time", y = "sigma") +
    .ts_theme()

  (p_mean / p_vol) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste0("Mean: ", sg_result$model_info$mean_model,
                        "  |  Variance: ", sg_result$model_info$variance_model,
                        " [", sg_result$model_info$garch_dist, "]"),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(colour = "grey40")
      )
    )
}


# =============================================================================
# 4. plot_ts_garch()  — Volatility Visualisation
# =============================================================================

#' GARCH Volatility Plot
#'
#' Four-panel overview: (1) observed returns with \eqn{\pm 2\sigma} bands, (2) conditional
#' variance, (3) standardised residuals, (4) histogram of standardised residuals.
#'
#' @param garch_result Result list from \code{ts_garch()}.
#' @param title Character. Plot title.
#' @importFrom rlang .data
#'
#' @return A \code{patchwork} composite \code{ggplot}.
#'
#' @examples
#' r = diff(log(AirPassengers))
#' res = ts_garch(r)
#' plot_ts_garch(res)
#'
#' @export
plot_ts_garch = function(garch_result,
                          title = "GARCH Volatility Analysis") {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required.")

  tbl = garch_result$fitted
  info = garch_result$model_info

  # Panel 1: observed + 2-sigma bands
  p1 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = -2 * .data$sigma, ymax = 2 * .data$sigma),
                         fill = .PALETTE$sigma, alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = .data$observed),
                       colour = .PALETTE$observed, linewidth = 0.7) +
    ggplot2::labs(subtitle = "Returns with +/- 2 sigma Band",
                  x = NULL, y = "Return") +
    .ts_theme()

  # Panel 2: conditional variance
  p2 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index, y = .data$variance)) +
    ggplot2::geom_area(fill = .PALETTE$sigma, alpha = 0.4) +
    ggplot2::geom_line(colour = .PALETTE$sigma, linewidth = 0.8) +
    ggplot2::labs(subtitle = "Conditional Variance h_t",
                  x = NULL, y = "Variance") +
    .ts_theme()

  # Panel 3: standardised residuals
  p3 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index, y = .data$std_residual)) +
    ggplot2::geom_line(colour = .PALETTE$fitted, linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = c(-3, 3),
                        colour = "red", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(subtitle = "Standardised Residuals z_t",
                  x = "Time", y = "z") +
    .ts_theme()

  # Panel 4: histogram of std residuals
  p4 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$std_residual)) +
    ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)),
                            bins = 30,
                            fill = .PALETTE$fitted, colour = "white",
                            alpha = 0.7) +
    ggplot2::stat_function(fun = stats::dnorm, colour = .PALETTE$forecast,
                           linewidth = 1) +
    ggplot2::labs(subtitle = "Density of z_t vs N(0,1)",
                  x = "z", y = "Density") +
    .ts_theme()

  (p1 + p2) / (p3 + p4) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = paste0(info$model_type, " [", info$distribution,
                        "]  |  AIC: ", info$aic),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                size = 14))
    )
}


# =============================================================================
# 5. plot_ts_residuals()  — Residual Diagnostic Plots
# =============================================================================

#' Residual Diagnostic Plots
#'
#' Four-panel diagnostic: (1) residuals over time, (2) ACF of residuals,
#' (3) Q-Q plot, (4) histogram vs normal density.
#'
#' @param model_result Result from any \code{ts_*()} function that contains a
#'   \code{$fitted} tibble with a \code{residual} column.
#' @param title Character. Plot title.
#' @param max_lag Integer. Maximum lag for ACF (default \code{30}).
#' @importFrom rlang .data
#'
#' @return A \code{patchwork} composite \code{ggplot}.
#'
#' @examples
#' data(AirPassengers)
#' fit = ts_sarima(log(AirPassengers))
#' plot_ts_residuals(fit)
#'
#' @export
plot_ts_residuals = function(model_result,
                              title   = "Residual Diagnostics",
                              max_lag = 30L) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required.")

  # Locate residual vector
  # Locate residual vector
  resid = if ("std_residual" %in% names(model_result$fitted)) {
    model_result$fitted$std_residual      # prefer standardised (GARCH)
  } else if ("residual" %in% names(model_result$fitted)) {
    model_result$fitted$residual
  } else {
    stop("No residual column found in model_result$fitted.")
  }

  n   = length(resid)
  idx = seq_len(n)

  # Panel 1: residuals over time
  p1 = ggplot2::ggplot(tibble::tibble(t = idx, r = resid),
                        ggplot2::aes(.data$t, .data$r)) +
    ggplot2::geom_line(colour = .PALETTE$observed, linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = 0,
                        colour = .PALETTE$ref_line, linetype = "dashed") +
    ggplot2::labs(subtitle = "Residuals vs Time", x = "Time", y = "Residual") +
    .ts_theme()

  # Panel 2: ACF via tidy computation
  acf_obj  = stats::acf(resid, lag.max = max_lag, plot = FALSE)
  ci_bound = stats::qnorm(0.975) / sqrt(n)
  acf_tbl  = tibble::tibble(lag = as.numeric(acf_obj$lag[-1]),
                             acf = as.numeric(acf_obj$acf[-1]))
  p2 = ggplot2::ggplot(acf_tbl, ggplot2::aes(x = .data$lag, y = .data$acf)) +
    ggplot2::geom_hline(yintercept = c(-ci_bound, ci_bound),
                        linetype = "dashed", colour = .PALETTE$forecast,
                        linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                          colour = .PALETTE$observed, linewidth = 0.8) +
    ggplot2::geom_point(colour = .PALETTE$observed, size = 1.5) +
    ggplot2::labs(subtitle = "ACF of Residuals",
                  x = "Lag", y = "ACF") +
    .ts_theme()

  # Panel 3: Q-Q plot
  qq_tbl = tibble::tibble(
    theoretical = stats::qnorm(stats::ppoints(n)),
    sample      = sort(resid)
  )
  p3 = ggplot2::ggplot(qq_tbl, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_point(colour = .PALETTE$fitted, alpha = 0.6, size = 1.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = .PALETTE$forecast, linetype = "dashed") +
    ggplot2::labs(subtitle = "Normal Q-Q Plot",
                  x = "Theoretical Quantiles",
                  y = "Sample Quantiles") +
    .ts_theme()

  # Panel 4: Histogram vs normal
  p4 = ggplot2::ggplot(tibble::tibble(r = resid), ggplot2::aes(.data$r)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 30,
                            fill = .PALETTE$fitted, colour = "white",
                            alpha = 0.7) +
    ggplot2::stat_function(
      fun  = stats::dnorm,
      args = list(mean = mean(resid, na.rm = TRUE),
                  sd   = stats::sd(resid, na.rm = TRUE)),
      colour = .PALETTE$forecast, linewidth = 1
    ) +
    ggplot2::labs(subtitle = "Histogram vs N(mu, sigma^2)",
                  x = "Residual", y = "Density") +
    .ts_theme()

  (p1 + p2) / (p3 + p4) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                size = 14))
    )
}


# =============================================================================
# 6. plot_ts_stl()  — STL Decomposition Visualisation
# =============================================================================

#' Plot STL Decomposition Components
#'
#' Stacked four-panel plot: observed, trend, seasonal, remainder.
#'
#' @param stl_result Result list from \code{ts_stl()}.
#' @param title Character. Plot title.
#' @importFrom rlang .data
#'
#' @return A \code{patchwork} composite \code{ggplot}.
#'
#' @examples
#' data(AirPassengers)
#' res = ts_stl(AirPassengers)
#' plot_ts_stl(res)
#'
#' @export
plot_ts_stl = function(stl_result,
                        title = "STL Decomposition") {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required.")

  tbl = stl_result$components
  str = stl_result$strength

  make_panel = function(y_col, subtitle, colour, ref = FALSE) {
    p = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$index, y = .data[[y_col]])) +
      ggplot2::geom_line(colour = colour, linewidth = 0.9) +
      ggplot2::labs(subtitle = subtitle, x = NULL, y = y_col) +
      .ts_theme()
    if (ref) p = p + ggplot2::geom_hline(yintercept = 0,
                                          colour = .PALETTE$ref_line,
                                          linetype = "dashed")
    p
  }

  p_obs = make_panel("observed", "Observed", .PALETTE$observed)
  p_tr  = make_panel("trend",    paste0("Trend  [strength = ",
                                         str$trend_strength, "]"),
                      .PALETTE$trend)
  p_sea = make_panel("seasonal", paste0("Seasonal  [strength = ",
                                         str$seasonal_strength, "]"),
                      .PALETTE$seasonal)
  p_rem = make_panel("remainder", "Remainder",
                      .PALETTE$remainder, ref = TRUE)

  # x-axis label only on bottom panel
  p_rem = p_rem + ggplot2::labs(x = "Time Index")

  (p_obs / p_tr / p_sea / p_rem) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                size = 14))
    )
}


# =============================================================================
# 7. plot_ts_acf() + plot_ts_pacf()  — ACF / PACF Plots
# =============================================================================

#' ACF and PACF Plots
#'
#' @param x A numeric vector, \code{ts} object, or tidy tibble/data.frame with
#'   a \code{value} column.
#' @param max_lag Integer. Maximum lag (default 40).
#' @param type Character. \code{"acf"} (default), \code{"pacf"}, or
#'   \code{"both"} (returns a patchwork panel).
#' @param title Character. Plot title.
#' @param diff Integer. Apply \code{diff()} this many times before plotting
#'   (default 0).
#' @importFrom rlang .data
#'
#' @return A \code{ggplot} or \code{patchwork} object.
#'
#' @examples
#' data(AirPassengers)
#' plot_ts_acf(log(AirPassengers), type = "both")
#'
#' @export
plot_ts_acf = function(x,
                        max_lag = 40L,
                        type    = c("acf", "pacf", "both"),
                        title   = NULL,
                        diff    = 0L) {
  type = match.arg(type)

  # Extract numeric vector
  x_vec = if (is.data.frame(x) && "value" %in% names(x)) x$value else as.numeric(x)
  for (i in seq_len(diff)) x_vec = base::diff(x_vec)
  n        = length(x_vec)
  ci_bound = stats::qnorm(0.975) / sqrt(n)

  .make_acf_panel = function(vec, lag_max, is_pacf) {
    fn  = if (is_pacf) stats::pacf else stats::acf
    obj = fn(vec, lag.max = lag_max, plot = FALSE)
    tbl = tibble::tibble(
      lag  = as.numeric(obj$lag[if (is_pacf) seq_along(obj$lag) else -1]),
      corr = as.numeric(obj$acf[if (is_pacf) seq_along(obj$acf) else -1])
    )
    sub = if (is_pacf) "Partial Autocorrelation (PACF)" else "Autocorrelation (ACF)"
    ggplot2::ggplot(tbl, ggplot2::aes(x = .data$lag, y = .data$corr)) +
      ggplot2::geom_hline(yintercept = c(-ci_bound, ci_bound),
                          linetype = "dashed",
                          colour = .PALETTE$forecast, linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                            colour = .PALETTE$observed, linewidth = 0.8) +
      ggplot2::geom_point(colour = .PALETTE$observed, size = 1.5) +
      ggplot2::labs(subtitle = sub, x = "Lag",
                    y = if (is_pacf) "PACF" else "ACF") +
      .ts_theme()
  }

  diff_label = if (diff > 0) paste0("  [diff=", diff, "]") else ""
  default_title = paste0(toupper(type), " Plot", diff_label)

  if (type == "acf") {
    p = .make_acf_panel(x_vec, max_lag, is_pacf = FALSE)
    return(p + ggplot2::ggtitle(if (!is.null(title)) title else default_title))
  }
  if (type == "pacf") {
    p = .make_acf_panel(x_vec, max_lag, is_pacf = TRUE)
    return(p + ggplot2::ggtitle(if (!is.null(title)) title else default_title))
  }

  # "both"
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' required for type='both'.")

  p_acf  = .make_acf_panel(x_vec, max_lag, is_pacf = FALSE)
  p_pacf = .make_acf_panel(x_vec, max_lag, is_pacf = TRUE)

  (p_acf / p_pacf) +
    patchwork::plot_annotation(
      title = if (!is.null(title)) title else paste0("ACF & PACF", diff_label),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                size = 14))
    )
}

#' PACF Plot (convenience alias)
#'
#' Shortcut for \code{plot_ts_acf(x, type = "pacf")}.
#'
#' @inheritParams plot_ts_acf
#' @export
plot_ts_pacf = function(x, max_lag = 40L, title = NULL, diff = 0L) {
  plot_ts_acf(x, max_lag = max_lag, type = "pacf", title = title, diff = diff)
}
