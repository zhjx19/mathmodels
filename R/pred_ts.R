# =============================================================================
# Lightweight time-series toolkit for mathmodels
# =============================================================================

# -----------------------------------------------------------------------------
# ts_df data structure
# -----------------------------------------------------------------------------

.ts_time_ok = function(x) {
  is.numeric(x) || inherits(x, "Date") || inherits(x, "POSIXt")
}

#' Lightweight Time-Series Data Frame
#'
#' Creates a simple single-series time-series data frame with fixed columns
#' `time` and `value`, plus `frequency` and `start` attributes.
#'
#' @param data A data frame.
#' @param time Name of the time/index column.
#' @param value Name of the numeric value column.
#' @param frequency Positive numeric scalar. Defaults to 1.
#' @param start Optional start metadata. If `NULL`, uses the first sorted time.
#'
#' @return A `ts_df` object, which is also a data frame.
#' @examples
#' df = data.frame(time = c(2, 1, 3), value = c(12, 10, 15))
#' ts_df(df)
#' @export
ts_df = function(data,
                 time = "time",
                 value = "value",
                 frequency = 1,
                 start = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  if (!time %in% names(data)) {
    stop("time column not found.", call. = FALSE)
  }
  if (!value %in% names(data)) {
    stop("value column not found.", call. = FALSE)
  }

  out = data[, c(time, value), drop = FALSE]
  names(out) = c("time", "value")

  if (!.ts_time_ok(out$time)) {
    stop("time must be numeric, Date, or POSIXct.", call. = FALSE)
  }
  if (!is.numeric(out$value)) {
    stop("value must be numeric.", call. = FALSE)
  }
  if (anyDuplicated(out$time)) {
    stop("time must not contain duplicate values.", call. = FALSE)
  }

  out = out[order(out$time), , drop = FALSE]
  row.names(out) = NULL

  if (is.null(start)) {
    start = out$time[1]
  }

  structure(
    out,
    frequency = frequency,
    start = start,
    class = c("ts_df", "data.frame")
  ) |>
    validate_ts_df()
}

#' Convert Common Inputs to `ts_df`
#'
#' Converts a `ts` object, numeric vector, explicit data-frame columns, or an
#' existing `ts_df` to the package's lightweight time-series structure.
#'
#' @param x A `ts`, numeric vector, data frame, or `ts_df`.
#' @param time,value Required column names when `x` is a data frame.
#' @param frequency Optional frequency metadata.
#' @param start Optional start metadata.
#'
#' @return A `ts_df` object.
#' @examples
#' as_ts_df(ts(1:4, frequency = 4))
#' as_ts_df(c(3, 5, 8))
#' as_ts_df(data.frame(month = 1:3, demand = c(10, 12, 15)), time = "month", value = "demand")
#' @export
as_ts_df = function(x,
                    time = NULL,
                    value = NULL,
                    frequency = NULL,
                    start = NULL) {
  if (inherits(x, "ts_df")) {
    if (!is.null(frequency)) attr(x, "frequency") = frequency
    if (!is.null(start)) attr(x, "start") = start
    return(validate_ts_df(x))
  }

  if (inherits(x, "ts") || (is.numeric(x) && !is.matrix(x))) {
    out = data.frame(time = seq_along(x), value = as.numeric(x))
    if (is.null(frequency)) {
      frequency = if (inherits(x, "ts")) stats::frequency(x) else 1
    }
    if (is.null(start)) {
      start = if (inherits(x, "ts")) stats::start(x) else 1
    }
    return(ts_df(out, frequency = frequency, start = start))
  }

  if (is.data.frame(x)) {
    if (is.null(time) || is.null(value)) {
      stop("time and value must be supplied for data.frame inputs.", call. = FALSE)
    }
    return(ts_df(
      x,
      time = time,
      value = value,
      frequency = if (is.null(frequency)) 1 else frequency,
      start = start
    ))
  }

  stop(
    "x must be a ts object, numeric vector, data.frame, or ts_df.",
    call. = FALSE
  )
}

#' Check Whether an Object Is a `ts_df`
#'
#' @param x Object to check.
#' @return `TRUE` or `FALSE`.
#' @examples
#' x = ts_df(data.frame(time = 1:3, value = c(10, 12, 15)))
#' is_ts_df(x)
#' @export
is_ts_df = function(x) {
  inherits(x, "ts_df") &&
    is.data.frame(x) &&
    identical(names(x), c("time", "value")) &&
    .ts_time_ok(x$time) &&
    is.numeric(x$value) &&
    !is.null(attr(x, "frequency")) &&
    !is.null(attr(x, "start"))
}

#' Validate a `ts_df`
#'
#' @param x A `ts_df` object.
#' @param require_complete If `TRUE`, first and last values must not be `NA`.
#' @param require_no_na If `TRUE`, no values may be `NA`.
#'
#' @return A sorted `ts_df` object.
#' @examples
#' x = ts_df(data.frame(time = 1:3, value = c(10, 12, 15)))
#' validate_ts_df(x)
#' @export
validate_ts_df = function(x,
                          require_complete = TRUE,
                          require_no_na = FALSE) {
  if (!inherits(x, "ts_df") || !is.data.frame(x)) {
    stop("x must be a ts_df.", call. = FALSE)
  }
  if (!identical(names(x), c("time", "value"))) {
    stop("ts_df must contain exactly time and value columns.", call. = FALSE)
  }
  if (!.ts_time_ok(x$time)) {
    stop("time must be numeric, Date, or POSIXct.", call. = FALSE)
  }
  if (!is.numeric(x$value)) {
    stop("value must be numeric.", call. = FALSE)
  }
  if (nrow(x) == 0) {
    stop("x must contain at least one row.", call. = FALSE)
  }
  if (anyDuplicated(x$time)) {
    stop("time must not contain duplicate values.", call. = FALSE)
  }

  freq = attr(x, "frequency")
  if (is.null(freq) || !is.numeric(freq) || length(freq) != 1 || freq <= 0) {
    stop("frequency must be a positive numeric scalar.", call. = FALSE)
  }

  start = attr(x, "start")
  if (is.null(start)) {
    stop("start attribute is required.", call. = FALSE)
  }

  if (require_complete && (is.na(x$value[1]) || is.na(x$value[nrow(x)]))) {
    stop("First and last values must not be NA.", call. = FALSE)
  }
  if (require_no_na && anyNA(x$value)) {
    stop("value must not contain NA.", call. = FALSE)
  }

  out = x[order(x$time), , drop = FALSE]
  row.names(out) = NULL
  structure(
    out,
    frequency = freq,
    start = start,
    class = c("ts_df", "data.frame")
  )
}

.ts_restore = function(x, time, value) {
  out = data.frame(time = time, value = value)
  structure(
    out,
    frequency = attr(x, "frequency"),
    start = attr(x, "start"),
    class = c("ts_df", "data.frame")
  ) |>
    validate_ts_df(require_complete = FALSE)
}

# -----------------------------------------------------------------------------
# Completion and interpolation
# -----------------------------------------------------------------------------

#' Complete Missing Time Points
#'
#' Adds missing time points on a regular sequence and leaves their values as
#' `NA`.
#'
#' @param x A `ts_df`.
#' @param by Step size passed to `seq()`.
#'
#' @return A completed `ts_df`.
#' @examples
#' x = ts_df(data.frame(time = c(1, 2, 4), value = c(10, 12, 18)))
#' complete_ts_df(x, by = 1)
#' @export
complete_ts_df = function(x, by) {
  x = validate_ts_df(x, require_complete = TRUE)
  if (missing(by)) {
    stop("by must be supplied.", call. = FALSE)
  }

  full_time = seq(min(x$time), max(x$time), by = by)
  idx = match(full_time, x$time)
  .ts_restore(x, full_time, x$value[idx])
}

#' Complete and Interpolate a `ts_df`
#'
#' Completes missing time points and fills internal missing values by linear or
#' spline interpolation. Endpoint values must already be observed.
#'
#' @param x A `ts_df`.
#' @param by Step size passed to `seq()`.
#' @param method Interpolation method: `"linear"` or `"spline"`.
#'
#' @return An interpolated `ts_df`.
#' @examples
#' x = ts_df(data.frame(time = c(1, 3, 5), value = c(10, 20, 30)))
#' impute_ts_df(x, by = 1, method = "linear")
#' @export
impute_ts_df = function(x,
                         by,
                         method = c("linear", "spline")) {
  method = match.arg(method)
  x = complete_ts_df(x, by = by)

  ok = !is.na(x$value)
  if (!ok[1] || !ok[length(ok)]) {
    stop("First and last values must not be NA for interpolation.", call. = FALSE)
  }

  pos = seq_along(x$value)
  value = switch(
    method,
    linear = stats::approx(pos[ok], x$value[ok], xout = pos)$y,
    spline = stats::spline(pos[ok], x$value[ok], xout = pos)$y
  )

  .ts_restore(x, x$time, value) |>
    validate_ts_df(require_no_na = TRUE)
}

#' Drop Missing Values from a `ts_df`
#'
#' @param x A `ts_df`.
#' @return A `ts_df` with rows whose `value` is `NA` removed.
#' @examples
#' x = ts_df(data.frame(time = 1:4, value = c(10, NA, 12, 13)))
#' drop_na_ts_df(x)
#' @export
drop_na_ts_df = function(x) {
  x = validate_ts_df(x, require_complete = FALSE)
  .ts_restore(x, x$time[!is.na(x$value)], x$value[!is.na(x$value)])
}

# -----------------------------------------------------------------------------
# Internal adapters
# -----------------------------------------------------------------------------

.ts_df_to_ts = function(x) {
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  stats::ts(x$value, frequency = attr(x, "frequency"), start = attr(x, "start"))
}

.ts_tidy_fitted = function(x, fitted, residual) {
  tibble::tibble(
    time = x$time,
    observed = x$value,
    fitted = as.numeric(fitted),
    residual = as.numeric(residual)
  )
}

# -----------------------------------------------------------------------------
# Decomposition and diagnostics
# -----------------------------------------------------------------------------

#' STL Decomposition for a `ts_df`
#'
#' Decomposes a complete seasonal `ts_df` with `stats::stl()`.
#'
#' @param x A complete seasonal `ts_df` with no missing values.
#' @param s_window Seasonal smoothing window.
#' @param robust Use robust STL fitting.
#' @param ... Additional arguments passed to `stats::stl()`.
#'
#' @return A list with `components`, `strength`, and `model`.
#' @examples
#' ts_stl(as_ts_df(AirPassengers))
#' @export
ts_stl = function(x,
                  s_window = "periodic",
                  robust = TRUE,
                  ...) {
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  x_ts = .ts_df_to_ts(x)
  if (stats::frequency(x_ts) <= 1) {
    stop("STL decomposition requires frequency > 1.", call. = FALSE)
  }

  fit = stats::stl(x_ts, s.window = s_window, robust = robust, ...)
  comp = fit$time.series

  components = tibble::tibble(
    time = x$time,
    observed = x$value,
    trend = as.numeric(comp[, "trend"]),
    seasonal = as.numeric(comp[, "seasonal"]),
    remainder = as.numeric(comp[, "remainder"])
  )

  var_rem = stats::var(components$remainder, na.rm = TRUE)
  seasonal_strength = max(
    0,
    1 - var_rem / stats::var(components$seasonal + components$remainder,
                             na.rm = TRUE)
  )
  trend_strength = max(
    0,
    1 - var_rem / stats::var(components$trend + components$remainder,
                             na.rm = TRUE)
  )

  list(
    components = components,
    strength = tibble::tibble(
      seasonal_strength = round(seasonal_strength, 4),
      trend_strength = round(trend_strength, 4)
    ),
    model = fit
  )
}

#' Stationarity Tests for a Time Series
#'
#' Runs ADF and KPSS tests on a complete `ts_df`.
#'
#' @param x A complete `ts_df` with no missing values.
#' @param adf_lags Optional lag for the ADF test.
#'
#' @return A tibble with test results.
#' @examples
#' ts_test(as_ts_df(log(AirPassengers)))
#' @export
ts_test = function(x, adf_lags = NULL) {
  x_ts = .ts_df_to_ts(x)
  n = length(x_ts)

  adf_lag = if (is.null(adf_lags)) trunc((n - 1)^(1 / 3)) else adf_lags
  adf_res = suppressWarnings(tseries::adf.test(x_ts, k = adf_lag))
  kpss_res = suppressWarnings(tseries::kpss.test(x_ts, null = "Level"))

  tibble::tibble(
    test = c("ADF", "KPSS"),
    null_hypothesis = c(
      "Unit root present (non-stationary)",
      "Series is stationary"
    ),
    statistic = c(unname(adf_res$statistic), unname(kpss_res$statistic)),
    p_value = c(adf_res$p.value, kpss_res$p.value),
    conclusion = c(
      ifelse(adf_res$p.value < 0.05, "Stationary", "Non-stationary"),
      ifelse(kpss_res$p.value >= 0.05, "Stationary", "Non-stationary")
    )
  )
}

# -----------------------------------------------------------------------------
# Transformation and back transformation
# -----------------------------------------------------------------------------

#' Transform a `ts_df`
#'
#' Applies a variance-stabilising transform and optional differencing. Returns a
#' transformed `ts_df` plus explicit parameters for back-transformation.
#'
#' @param x A complete `ts_df` with no missing values.
#' @param method One of `"none"`, `"log"`, or `"boxcox"`.
#' @param lambda Optional Box-Cox lambda.
#' @param diff Number of regular differences.
#' @param seasonal_diff Number of seasonal differences.
#'
#' @return A list with `transformed`, `params`, and `summary`.
#' @examples
#' ts_transform(as_ts_df(AirPassengers), method = "log", diff = 1)
#' @export
ts_transform = function(x,
                         method = c("none", "log", "boxcox"),
                         lambda = NULL,
                         diff = 0L,
                         seasonal_diff = 0L) {
  method = match.arg(method)
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  freq = attr(x, "frequency")

  values = x$value
  steps = character(0)
  lambda_used = NA_real_

  if (method == "log") {
    if (any(values <= 0)) {
      stop("Log transform requires all values > 0.", call. = FALSE)
    }
    values = log(values)
    steps = c(steps, "log")
  } else if (method == "boxcox") {
    lambda_used = if (!is.null(lambda)) lambda else forecast::BoxCox.lambda(values)
    values = forecast::BoxCox(values, lambda_used)
    steps = c(steps, "boxcox")
  }

  base_values = values
  drop_n = 0L

  if (seasonal_diff > 0L) {
    for (i in seq_len(seasonal_diff)) {
      values = diff(values, lag = freq)
      drop_n = drop_n + freq
    }
    steps = c(steps, "seasonal_diff")
  }

  if (diff > 0L) {
    for (i in seq_len(diff)) {
      values = diff(values)
      drop_n = drop_n + 1L
    }
    steps = c(steps, "diff")
  }

  transformed_time = x$time[(drop_n + 1L):nrow(x)]
  transformed = .ts_restore(x, transformed_time, as.numeric(values))

  params = list(
    method = method,
    lambda = lambda_used,
    diff = diff,
    seasonal_diff = seasonal_diff,
    frequency = freq,
    last_base_values = utils::tail(base_values, max(1L, freq))
  )

  summary = tibble::tibble(
    original_length = nrow(x),
    transformed_length = nrow(transformed),
    steps_applied = if (length(steps) == 0) "none" else paste(steps, collapse = " -> "),
    lambda = lambda_used
  )

  list(transformed = transformed, params = params, summary = summary)
}

#' Back-Transform Forecasts
#'
#' Converts forecasts produced on a transformed scale back to the original scale.
#'
#' @param forecast_tbl Forecast tibble from `ts_forecast()`.
#' @param params The `params` element from `ts_transform()`.
#'
#' @return The forecast tibble with numeric forecast columns restored.
#' @examples
#' x = as_ts_df(AirPassengers)
#' tr = ts_transform(x, method = "log")
#' fc = tibble::tibble(step = 1:2, forecast = log(c(500, 600)))
#' ts_back_transform(fc, tr$params)
#' @export
ts_back_transform = function(forecast_tbl, params) {
  tbl = forecast_tbl
  num_cols = setdiff(names(tbl)[vapply(tbl, is.numeric, logical(1))], "step")

  if (params$diff > 0L) {
    last_val = utils::tail(params$last_base_values, 1)
    for (col in num_cols) {
      tbl[[col]] = cumsum(c(last_val, tbl[[col]]))[-1]
    }
  }

  if (params$seasonal_diff > 0L) {
    freq = params$frequency
    last_season = utils::tail(params$last_base_values, freq)
    for (col in num_cols) {
      out = numeric(nrow(tbl))
      for (i in seq_len(nrow(tbl))) {
        base_val = if (i <= freq) last_season[i] else out[i - freq]
        out[i] = tbl[[col]][i] + base_val
      }
      tbl[[col]] = out
    }
  }

  if (params$method == "log") {
    for (col in num_cols) tbl[[col]] = exp(tbl[[col]])
  } else if (params$method == "boxcox") {
    for (col in num_cols) tbl[[col]] = forecast::InvBoxCox(tbl[[col]], params$lambda)
  }

  tbl
}

# -----------------------------------------------------------------------------
# Models and forecasts
# -----------------------------------------------------------------------------

#' ETS Model Fitting
#'
#' @param x A complete `ts_df` with no missing values.
#' @param model ETS model string passed to `forecast::ets()`.
#' @param ... Additional arguments passed to `forecast::ets()`.
#'
#' @return A list containing tidy model summaries, the raw model, and input.
#' @examples
#' ts_ets(as_ts_df(log(AirPassengers)))
#' @export
ts_ets = function(x, model = "ZZZ", ...) {
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  x_ts = .ts_df_to_ts(x)
  fit = forecast::ets(x_ts, model = model, ...)

  pars = c(fit$par)
  list(
    model_info = tibble::tibble(
      model_type = fit$method,
      log_lik = round(fit$loglik, 4),
      aic = round(fit$aic, 4),
      aicc = round(fit$aicc, 4),
      bic = round(fit$bic, 4)
    ),
    parameters = tibble::tibble(
      term = names(pars),
      estimate = round(as.numeric(pars), 6)
    ),
    fitted = .ts_tidy_fitted(x, fit$fitted, fit$residuals),
    model = fit,
    input = x
  )
}

#' SARIMA Model Fitting
#'
#' @param x A complete `ts_df` with no missing values.
#' @param order Optional non-seasonal ARIMA order.
#' @param seasonal Optional seasonal specification.
#' @param stepwise,approximation Passed to `forecast::auto.arima()`.
#' @param ... Additional arguments passed to the forecast fitting function.
#'
#' @return A list containing tidy model summaries, the raw model, and input.
#' @examples
#' ts_sarima(as_ts_df(log(AirPassengers)), stepwise = TRUE, approximation = TRUE)
#' @export
ts_sarima = function(x,
                      order = NULL,
                      seasonal = NULL,
                      stepwise = TRUE,
                      approximation = TRUE,
                      ...) {
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  x_ts = .ts_df_to_ts(x)

  if (is.null(order)) {
    fit = forecast::auto.arima(
      x_ts,
      stepwise = stepwise,
      approximation = approximation,
      ...
    )
  } else {
    seas = if (is.null(seasonal)) list(order = c(0, 0, 0)) else seasonal
    fit = forecast::Arima(x_ts, order = order, seasonal = seas, ...)
  }

  coef_vec = stats::coef(fit)
  se_vec = if (length(coef_vec) > 0 && !is.null(fit$var.coef)) {
    sqrt(diag(fit$var.coef))
  } else {
    numeric(0)
  }

  lb = stats::Box.test(
    fit$residuals,
    lag = min(20, floor(length(x_ts) / 5)),
    type = "Ljung-Box",
    fitdf = length(coef_vec)
  )

  list(
    model_info = tibble::tibble(
      model_type = fit$method,
      log_lik = round(fit$loglik, 4),
      aic = round(fit$aic, 4),
      aicc = round(fit$aicc, 4),
      bic = round(fit$bic, 4),
      sigma2 = round(fit$sigma2, 6)
    ),
    coefficients = tibble::tibble(
      term = names(coef_vec),
      estimate = round(as.numeric(coef_vec), 6),
      std_error = round(as.numeric(se_vec), 6)
    ),
    fitted = .ts_tidy_fitted(x, fit$fitted, fit$residuals),
    diagnostics = tibble::tibble(
      test = "Ljung-Box",
      lag = unname(lb$parameter),
      statistic = round(unname(lb$statistic), 4),
      p_value = round(lb$p.value, 4),
      conclusion = ifelse(lb$p.value > 0.05,
                          "No autocorrelation detected",
                          "Autocorrelation detected")
    ),
    model = fit,
    input = x
  )
}

#' Generate Forecasts
#'
#' @param model_result Result from `ts_ets()` or `ts_sarima()`.
#' @param h Forecast horizon.
#' @param level Confidence levels.
#'
#' @return A tibble with point forecasts and intervals.
#' @examples
#' fit = ts_ets(as_ts_df(log(AirPassengers)))
#' ts_forecast(fit, h = 3)
#' @export
ts_forecast = function(model_result, h = 12, level = c(80, 95)) {
  if (!inherits(model_result$model, c("ets", "ARIMA", "Arima"))) {
    stop("model_result must be returned by ts_ets() or ts_sarima().", call. = FALSE)
  }

  fc = forecast::forecast(model_result$model, h = h, level = level)
  tbl = tibble::tibble(step = seq_len(h), forecast = as.numeric(fc$mean))

  for (lv in level) {
    idx = which(fc$level == lv)
    tbl[[paste0("lo_", lv)]] = as.numeric(fc$lower[, idx])
    tbl[[paste0("hi_", lv)]] = as.numeric(fc$upper[, idx])
  }

  tbl
}

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

.ts_theme = function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

.acf_like_tbl = function(x, max_lag = 40, type = c("acf", "pacf")) {
  type = match.arg(type)
  x = validate_ts_df(x, require_complete = TRUE, require_no_na = TRUE)
  vec = x$value
  n = length(vec)
  ci = stats::qnorm(0.975) / sqrt(n)

  obj = if (type == "acf") {
    stats::acf(vec, lag.max = max_lag, plot = FALSE)
  } else {
    stats::pacf(vec, lag.max = max_lag, plot = FALSE)
  }

  lags = as.numeric(obj$lag)
  vals = as.numeric(obj$acf)
  if (type == "acf") {
    lags = lags[-1]
    vals = vals[-1]
  }

  tibble::tibble(lag = lags, value = vals, ci = ci)
}

.acf_like_plot = function(tbl, title, y_lab) {
  ggplot2::ggplot(tbl, ggplot2::aes(x = .data$lag, y = .data$value)) +
    ggplot2::geom_hline(yintercept = c(-tbl$ci[1], tbl$ci[1]),
                        linetype = "dashed", colour = "#E74C3C", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                          colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#2C3E50", size = 1.5) +
    ggplot2::labs(title = title, x = "Lag", y = y_lab) +
    .ts_theme()
}

.future_time = function(time, h, by) {
  if (missing(by) || is.null(by)) {
    if (is.numeric(time)) return(max(time) + seq_len(h))
    stop("by must be supplied when time is Date or POSIXct.", call. = FALSE)
  }
  seq(max(time), by = by, length.out = h + 1L)[-1L]
}

#' Plot a Time Series
#'
#' @param x A `ts_df`.
#' @param title Plot title.
#' @param x_lab,y_lab Axis labels.
#'
#' @return A ggplot object.
#' @examples
#' plot_ts(as_ts_df(AirPassengers))
#' @export
plot_ts = function(x,
                    title = "Time Series",
                    x_lab = "Time",
                    y_lab = "Value") {
  x = validate_ts_df(x, require_complete = FALSE)

  ggplot2::ggplot(x, ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    .ts_theme()
}

#' Autocorrelation Plot
#'
#' @param x A complete `ts_df` with no missing values.
#' @param max_lag Maximum lag.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @examples
#' plot_ts_acf(as_ts_df(log(AirPassengers)))
#' @export
plot_ts_acf = function(x, max_lag = 40L, title = "ACF") {
  tbl = .acf_like_tbl(x, max_lag = max_lag, type = "acf")
  .acf_like_plot(tbl, title = title, y_lab = "ACF")
}

#' Partial Autocorrelation Plot
#'
#' @param x A complete `ts_df` with no missing values.
#' @param max_lag Maximum lag.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @examples
#' plot_ts_pacf(as_ts_df(log(AirPassengers)))
#' @export
plot_ts_pacf = function(x, max_lag = 40L, title = "PACF") {
  tbl = .acf_like_tbl(x, max_lag = max_lag, type = "pacf")
  .acf_like_plot(tbl, title = title, y_lab = "PACF")
}

#' STL Decomposition Plot
#'
#' @param stl_result Result from `ts_stl()`.
#' @param title Plot title.
#'
#' @return A patchwork plot.
#' @examples
#' plot_ts_stl(ts_stl(as_ts_df(AirPassengers)))
#' @export
plot_ts_stl = function(stl_result, title = "STL Decomposition") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork').", call. = FALSE)
  }

  tbl = stl_result$components
  str = stl_result$strength

  make_panel = function(y_col, subtitle, colour, ref = FALSE) {
    p = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$time, y = .data[[y_col]])) +
      ggplot2::geom_line(colour = colour, linewidth = 0.9) +
      ggplot2::labs(subtitle = subtitle, x = NULL, y = y_col) +
      .ts_theme()
    if (ref) {
      p = p + ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed")
    }
    p
  }

  p_obs = make_panel("observed", "Observed", "#2C3E50")
  p_tr  = make_panel("trend", paste0("Trend [", str$trend_strength, "]"), "#27AE60")
  p_sea = make_panel("seasonal", paste0("Seasonal [", str$seasonal_strength, "]"), "#F39C12")
  p_rem = make_panel("remainder", "Remainder", "#95A5A6", ref = TRUE) +
    ggplot2::labs(x = "Time")

  (p_obs / p_tr / p_sea / p_rem) +
    patchwork::plot_annotation(title = title)
}

#' Residual Diagnostic Plot
#'
#' @param model_result Result from `ts_ets()` or `ts_sarima()`.
#' @param max_lag Maximum lag.
#' @param title Plot title.
#'
#' @return A patchwork plot.
#' @examples
#' fit = ts_sarima(as_ts_df(log(AirPassengers)), stepwise = TRUE, approximation = TRUE)
#' plot_ts_residuals(fit)
#' @export
plot_ts_residuals = function(model_result,
                              max_lag = 30L,
                              title = "Residual Diagnostics") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork').", call. = FALSE)
  }

  resid = model_result$fitted$residual
  tbl = tibble::tibble(time = model_result$fitted$time, residual = resid)

  p1 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$time, y = .data$residual)) +
    ggplot2::geom_line(colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    ggplot2::labs(subtitle = "Residuals over time", x = NULL, y = "Residual") +
    .ts_theme()

  acf_tbl = .acf_like_tbl(structure(data.frame(time = tbl$time, value = tbl$residual),
                                    frequency = 1, start = 1,
                                    class = c("ts_df", "data.frame")),
                          max_lag = max_lag, type = "acf")
  p2 = .acf_like_plot(acf_tbl, title = NULL, y_lab = "ACF") +
    ggplot2::labs(subtitle = "Residual ACF")

  p3 = ggplot2::ggplot(tbl, ggplot2::aes(sample = .data$residual)) +
    ggplot2::stat_qq(colour = "#2C3E50") +
    ggplot2::stat_qq_line(colour = "#E74C3C") +
    ggplot2::labs(subtitle = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    .ts_theme()

  p4 = ggplot2::ggplot(tbl, ggplot2::aes(x = .data$residual)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat("density")),
                            bins = 30, fill = "#2C3E50", colour = "white", alpha = 0.7) +
    ggplot2::stat_function(fun = stats::dnorm,
                           args = list(mean = mean(resid, na.rm = TRUE), sd = stats::sd(resid, na.rm = TRUE)),
                           colour = "#E74C3C") +
    ggplot2::labs(subtitle = "Histogram", x = "Residual", y = "Density") +
    .ts_theme()

  (p1 + p2) / (p3 + p4) + patchwork::plot_annotation(title = title)
}

#' Plot Forecasts
#'
#' @param x Original `ts_df`.
#' @param forecast_tbl Result from `ts_forecast()`.
#' @param level Confidence levels to draw when present.
#' @param by Step size for future time values.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @examples
#' fit = ts_ets(as_ts_df(log(AirPassengers)))
#' fc = ts_forecast(fit, h = 3)
#' plot_ts_forecast(as_ts_df(log(AirPassengers)), fc, by = 1)
#' @export
plot_ts_forecast = function(x,
                             forecast_tbl,
                             level = c(80, 95),
                             by = NULL,
                             title = "Forecast") {
  x = validate_ts_df(x, require_complete = FALSE)
  future = .future_time(x$time, nrow(forecast_tbl), by = by)
  fc = forecast_tbl
  fc$time = future

  p = plot_ts(x, title = NULL)

  for (lv in sort(level, decreasing = TRUE)) {
    lo = paste0("lo_", lv)
    hi = paste0("hi_", lv)
    if (lo %in% names(fc) && hi %in% names(fc)) {
      p = p + ggplot2::geom_ribbon(
        data = fc,
        ggplot2::aes(x = .data$time, ymin = .data[[lo]], ymax = .data[[hi]]),
        fill = "#E74C3C",
        alpha = ifelse(lv == max(level), 0.12, 0.20),
        inherit.aes = FALSE
      )
    }
  }

  p +
    ggplot2::geom_line(
      data = fc,
      ggplot2::aes(x = .data$time, y = .data$forecast),
      colour = "#E74C3C",
      linewidth = 0.9,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::labs(title = title, x = "Time", y = "Value") +
    .ts_theme()
}
