# =============================================================================
# Unit tests for ts_models.R and ts_plots.R
# =============================================================================

# ---- Helpers ----
data(AirPassengers)
ap_ts = AirPassengers
ap_log = log(ap_ts)
ap_ret = diff(ap_log)

# =============================================================================
# ts_test() — Stationarity Tests
# =============================================================================

test_that("ts_test works with ts object", {
  res = ts_test(ap_log, verbose = FALSE)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 3)
  expect_setequal(res$test, c("ADF", "KPSS", "PP"))
  expect_type(res$statistic, "double")
  expect_type(res$p_value, "double")
  expect_true(all(res$conclusion %in% c("Stationary", "Non-stationary")))
})

test_that("ts_test works with numeric vector", {
  res = ts_test(as.numeric(ap_log), verbose = FALSE)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 3)
})

test_that("ts_test works with data.frame", {
  df = data.frame(val = as.numeric(ap_log))
  res = ts_test(df, verbose = FALSE)
  expect_s3_class(res, "tbl_df")
})

test_that("ts_test verbose=TRUE prints output", {
  expect_output(ts_test(ap_log, verbose = TRUE), "Stationarity Tests")
})

test_that("ts_test KPSS type 'Trend' works", {
  res = ts_test(ap_log, kpss_type = "Trend", verbose = FALSE)
  expect_equal(res$test[2], "KPSS")
})

test_that("ts_test PP type 'Z(alpha)' works", {
  res = ts_test(ap_log, pp_type = "Z(alpha)", verbose = FALSE)
  expect_equal(res$test[3], "PP")
})

# =============================================================================
# ts_stl() — STL Decomposition
# =============================================================================
# =============================================================================
# ts_stl() — STL Decomposition
# =============================================================================

test_that("ts_stl works with AirPassengers", {
  res = ts_stl(ap_ts)
  expect_type(res, "list")
  expect_setequal(names(res), c("components", "strength", "model"))
  expect_s3_class(res$components, "tbl_df")
  expect_equal(nrow(res$components), length(ap_ts))
  expect_setequal(names(res$components),
                  c("index", "observed", "trend", "seasonal", "remainder"))
  expect_s3_class(res$strength, "tbl_df")
  expect_setequal(names(res$strength),
                  c("seasonal_strength", "trend_strength"))
  expect_true(res$strength$seasonal_strength >= 0)
  expect_true(res$strength$trend_strength >= 0)
})

test_that("ts_stl works with numeric + frequency", {
  x_num = as.numeric(ap_ts)
  res = ts_stl(x_num, frequency = 12L)
  expect_equal(nrow(res$components), length(x_num))
})

test_that("ts_stl errors with frequency <= 1", {
  x = ts(1:20, frequency = 1)
  expect_error(ts_stl(x), "frequency > 1")
})

test_that("ts_stl s_window = 'periodic' by default", {
  res = ts_stl(ap_ts)
  expect_s3_class(res$model, "stl")
})

test_that("ts_stl robust = FALSE works", {
  res = ts_stl(ap_ts, robust = FALSE)
  expect_s3_class(res$components, "tbl_df")
})

# =============================================================================
# ts_ets() — Exponential Smoothing
# =============================================================================

test_that("ts_ets works with auto model (ZZZ)", {
  res = ts_ets(ap_log)
  expect_type(res, "list")
  expect_setequal(names(res), c("model_info", "parameters", "fitted", "model"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_equal(nrow(res$model_info), 1)
  expect_setequal(names(res$model_info),
                  c("model_type", "log_lik", "aic", "aicc", "bic"))
  expect_type(res$model_info$aic, "double")
})

test_that("ts_ets returns fitted tibble with correct columns", {
  res = ts_ets(ap_log)
  expect_s3_class(res$fitted, "tbl_df")
  expect_equal(nrow(res$fitted), length(ap_log))
  expect_setequal(names(res$fitted),
                  c("index", "observed", "fitted", "residual"))
})

test_that("ts_ets works with user-specified model 'AAN'", {
  res = ts_ets(ap_log, model = "AAN")
  expect_match(res$model_info$model_type, "ETS\\(A,A,N\\)")
})

# =============================================================================
# ts_sarima() — SARIMA Modeling
# =============================================================================

test_that("ts_sarima works with auto selection", {
  res = ts_sarima(ap_log)
  expect_type(res, "list")
  expect_setequal(names(res),
                  c("model_info", "coefficients", "fitted", "diagnostics", "model"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_equal(nrow(res$model_info), 1)
  expect_setequal(names(res$model_info),
                  c("model_type", "log_lik", "aic", "aicc", "bic", "sigma2"))
})

test_that("ts_sarima auto-arima model_type starts with 'ARIMA'", {
  res = ts_sarima(ap_log)
  expect_match(res$model_info$model_type, "^ARIMA")
})

test_that("ts_sarima coefficients have expected columns", {
  res = ts_sarima(ap_log)
  expect_setequal(names(res$coefficients),
                  c("term", "estimate", "std_error", "t_stat"))
})

test_that("ts_sarima fitted tibble is correct length", {
  res = ts_sarima(ap_log)
  expect_equal(nrow(res$fitted), length(ap_log))
  expect_setequal(names(res$fitted),
                  c("index", "observed", "fitted", "residual"))
})

test_that("ts_sarima diagnostics returns Ljung-Box", {
  res = ts_sarima(ap_log)
  expect_s3_class(res$diagnostics, "tbl_df")
  expect_equal(nrow(res$diagnostics), 1)
  expect_equal(res$diagnostics$test, "Ljung-Box")
  expect_true(res$diagnostics$conclusion %in%
              c("No autocorrelation (residuals OK)", "Autocorrelation detected"))
})

test_that("ts_sarima works with user-specified order", {
  res = ts_sarima(ap_log, order = c(2, 1, 1))
  expect_match(res$model_info$model_type, "ARIMA2,1,1")
})

# =============================================================================
# ts_garch() — GARCH Variance Modeling
# =============================================================================

test_that("ts_garch works with return series", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret)
  expect_type(res, "list")
  expect_setequal(names(res),
                  c("model_info", "coefficients", "fitted", "diagnostics", "model"))
})

test_that("ts_garch model_info has expected columns", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret)
  expect_s3_class(res$model_info, "tbl_df")
  expect_equal(nrow(res$model_info), 1)
  expect_setequal(names(res$model_info),
                  c("model_type", "distribution", "log_lik", "aic", "bic"))
})

test_that("ts_garch fitted tibble has sigma and variance", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret)
  expect_s3_class(res$fitted, "tbl_df")
  expect_equal(nrow(res$fitted), length(ap_ret))
  expect_setequal(names(res$fitted),
                  c("index", "observed", "sigma", "variance",
                    "residual", "std_residual"))
})

test_that("ts_garch diagnostics returns ARCH-LM and Ljung-Box", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret)
  expect_s3_class(res$diagnostics, "tbl_df")
  expect_equal(nrow(res$diagnostics), 2)
  expect_true(any(grepl("ARCH-LM", res$diagnostics$test)))
  expect_true(any(grepl("Ljung-Box", res$diagnostics$test)))
})

test_that("ts_garch works with Student-t distribution", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret, dist = "std")
  expect_equal(res$model_info$distribution, "std")
})

test_that("ts_garch coefficients have expected columns", {
  skip_if_not_installed("rugarch")
  res = ts_garch(ap_ret)
  expect_setequal(names(res$coefficients),
                  c("term", "estimate", "std_error", "t_stat", "p_value"))
})

# =============================================================================
# ts_sarima_garch() — Two-stage SARIMA+GARCH
# =============================================================================

test_that("ts_sarima_garch returns combined result", {
  skip_if_not_installed("rugarch")
  res = ts_sarima_garch(ap_log)
  expect_type(res, "list")
  expect_setequal(names(res),
                  c("mean_model", "variance_model", "fitted", "model_info"))
})

test_that("ts_sarima_garch combined fitted has expected columns", {
  skip_if_not_installed("rugarch")
  res = ts_sarima_garch(ap_log)
  expect_s3_class(res$fitted, "tbl_df")
  expect_equal(nrow(res$fitted), length(ap_log))
  expect_setequal(names(res$fitted),
                  c("index", "observed", "mean_fitted",
                    "sigma", "variance", "std_residual"))
})

test_that("ts_sarima_garch model_info combines both models", {
  skip_if_not_installed("rugarch")
  res = ts_sarima_garch(ap_log)
  expect_s3_class(res$model_info, "tbl_df")
  expect_equal(nrow(res$model_info), 1)
  expect_setequal(names(res$model_info),
                  c("mean_model", "variance_model", "garch_dist",
                    "sarima_aic", "garch_aic"))
})

# =============================================================================
# ts_forecast() — Generic Forecast
# =============================================================================

test_that("ts_forecast works with SARIMA result", {
  fit = ts_sarima(ap_log)
  fc = ts_forecast(fit, h = 6)
  expect_s3_class(fc, "tbl_df")
  expect_equal(nrow(fc), 6)
  expect_setequal(names(fc), c("step", "forecast", "lo_80", "hi_80", "lo_95", "hi_95"))
  expect_type(fc$forecast, "double")
})

test_that("ts_forecast works with ETS result", {
  fit = ts_ets(ap_log)
  fc = ts_forecast(fit, h = 6)
  expect_equal(nrow(fc), 6)
  expect_true(all(c("lo_80", "hi_80") %in% names(fc)))
})

test_that("ts_forecast works with GARCH result", {
  skip_if_not_installed("rugarch")
  fit = ts_garch(ap_ret)
  fc = ts_forecast(fit, h = 6)
  expect_equal(nrow(fc), 6)
  expect_true("sigma" %in% names(fc))
})

test_that("ts_forecast works with SARIMA-GARCH result", {
  skip_if_not_installed("rugarch")
  fit = ts_sarima_garch(ap_log)
  fc = ts_forecast(fit, h = 6)
  expect_equal(nrow(fc), 6)
  expect_true("sigma" %in% names(fc))
})

test_that("ts_forecast custom levels work", {
  fit = ts_sarima(ap_log)
  fc = ts_forecast(fit, h = 6, level = c(90, 99))
  expect_true(all(c("lo_90", "hi_90", "lo_99", "hi_99") %in% names(fc)))
})

test_that("ts_forecast errors on unrecognised model", {
  expect_error(ts_forecast(list(model = "not_a_model"), h = 6),
               "Unrecognised")
})

# =============================================================================
# plot_ts() — Basic Time Series Plot
# =============================================================================

test_that("plot_ts works with ts object", {
  p = plot_ts(ap_ts, title = "Test")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts works with numeric vector", {
  p = plot_ts(as.numeric(ap_ts))
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts works with tidy tibble", {
  tbl = tibble::tibble(index = 1:30, value = rnorm(30))
  p = plot_ts(tbl)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts works with named list (multi-series)", {
  skip_if_not_installed("purrr")
  x_list = list(series_a = rnorm(30), series_b = rnorm(30))
  p = plot_ts(x_list)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts add_points = TRUE adds points", {
  p1 = plot_ts(ap_ts, add_points = FALSE)
  p2 = plot_ts(ap_ts, add_points = TRUE)
  expect_equal(length(p1$layers), 1)  # line only
  expect_equal(length(p2$layers), 2)  # line + point
})

test_that("plot_ts custom labels work", {
  p = plot_ts(ap_ts, title = "MyTitle", subtitle = "MySub",
              x_lab = "X", y_lab = "Y")
  expect_equal(p$labels$title, "MyTitle")
  expect_equal(p$labels$subtitle, "MySub")
})

# =============================================================================
# plot_ts_forecast() — Forecast Plot
# =============================================================================

test_that("plot_ts_forecast works with SARIMA forecast", {
  fit = ts_sarima(ap_log)
  fc  = ts_forecast(fit, h = 24)
  p   = plot_ts_forecast(ap_log, fc, title = "Test")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts_forecast trim_n limits history", {
  fit = ts_sarima(ap_log)
  fc  = ts_forecast(fit, h = 12)
  p   = plot_ts_forecast(ap_log, fc, trim_n = 24)
  expect_s3_class(p, "ggplot")
})

# =============================================================================
# plot_ts_sarima_garch() — SARIMA+GARCH Dual Panel
# =============================================================================

test_that("plot_ts_sarima_garch returns patchwork", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("patchwork")
  res = ts_sarima_garch(ap_log)
  p = plot_ts_sarima_garch(res)
  expect_s3_class(p, "patchwork")
})

# =============================================================================
# plot_ts_garch() — GARCH Volatility Panel
# =============================================================================

test_that("plot_ts_garch returns patchwork", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("patchwork")
  res = ts_garch(ap_ret)
  p = plot_ts_garch(res)
  expect_s3_class(p, "patchwork")
})

# =============================================================================
# plot_ts_residuals() — Residual Diagnostics
# =============================================================================

test_that("plot_ts_residuals works with SARIMA result", {
  skip_if_not_installed("patchwork")
  fit = ts_sarima(ap_log)
  p = plot_ts_residuals(fit)
  expect_s3_class(p, "patchwork")
})

test_that("plot_ts_residuals prefers std_residual from GARCH", {
  skip_if_not_installed("rugarch")
  skip_if_not_installed("patchwork")
  res = ts_garch(ap_ret)
  p = plot_ts_residuals(res)
  expect_s3_class(p, "patchwork")
})

# =============================================================================
# plot_ts_stl() — STL Decomposition Plot
# =============================================================================

test_that("plot_ts_stl returns patchwork", {
  skip_if_not_installed("patchwork")
  res = ts_stl(ap_ts)
  p = plot_ts_stl(res)
  expect_s3_class(p, "patchwork")
})

# =============================================================================
# plot_ts_acf() / plot_ts_pacf() — ACF/PACF
# =============================================================================

test_that("plot_ts_acf type='acf' works", {
  p = plot_ts_acf(ap_ret, type = "acf")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts_acf type='pacf' works", {
  p = plot_ts_acf(ap_ret, type = "pacf")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts_acf type='both' returns patchwork", {
  skip_if_not_installed("patchwork")
  p = plot_ts_acf(ap_ret, type = "both")
  expect_s3_class(p, "patchwork")
})

test_that("plot_ts_acf with diff works", {
  p = plot_ts_acf(ap_log, type = "acf", diff = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ts_pacf convenience alias works", {
  p1 = plot_ts_pacf(ap_ret)
  p2 = plot_ts_acf(ap_ret, type = "pacf")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_ts_acf works with data.frame input", {
  df = data.frame(value = as.numeric(ap_ret))
  p = plot_ts_acf(df, type = "acf")
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# ts_transform() — Preprocessing: Box-Cox, Log, Differencing
# =============================================================================

test_that("ts_transform method='none' returns identity", {
  res = ts_transform(ap_ts, method = "none", diff = 0, seasonal_diff = 0,
                     verbose = FALSE)
  expect_equal(as.numeric(res$transformed), as.numeric(ap_ts))
  expect_equal(res$params$method, "none")
  expect_equal(res$params$diff, 0L)
  expect_equal(res$params$seasonal_diff, 0L)
  expect_equal(res$summary$steps_applied, "none")
})

test_that("ts_transform method='log' works", {
  res = ts_transform(ap_ts, method = "log", diff = 0, seasonal_diff = 0,
                     verbose = FALSE)
  expect_equal(as.numeric(res$transformed), as.numeric(log(ap_ts)))
  expect_equal(res$params$method, "log")
  expect_equal(res$summary$steps_applied, "log(x)")
})

test_that("ts_transform method='boxcox' auto-selects lambda", {
  skip_if(!requireNamespace("forecast", quietly = TRUE),
          "forecast required")
  res = ts_transform(ap_ts, method = "boxcox", diff = 0, seasonal_diff = 0,
                     verbose = FALSE)
  expect_true(!is.na(res$params$lambda))
  expect_match(res$summary$steps_applied, "BoxCox")
})

test_that("ts_transform method='boxcox' with fixed lambda", {
  skip_if(!requireNamespace("forecast", quietly = TRUE),
          "forecast required")
  res = ts_transform(ap_ts, method = "boxcox", lambda = 0, diff = 0,
                     seasonal_diff = 0, verbose = FALSE)
  expect_equal(res$params$lambda, 0)
})

test_that("ts_transform log errors on non-positive values", {
  x = c(-1, 1, 2, 3)
  expect_error(ts_transform(x, method = "log", verbose = FALSE),
               "Log transform requires")
})

test_that("ts_transform diff=1 works", {
  res = ts_transform(ap_ts, method = "none", diff = 1, seasonal_diff = 0,
                     verbose = FALSE)
  expect_equal(length(res$transformed), length(ap_ts) - 1)
  expect_equal(res$params$diff, 1L)
})

test_that("ts_transform diff=2 works", {
  res = ts_transform(ap_ts, method = "none", diff = 2, seasonal_diff = 0,
                     verbose = FALSE)
  expect_equal(length(res$transformed), length(ap_ts) - 2)
  expect_equal(res$params$diff, 2L)
})

test_that("ts_transform seasonal_diff=1 works", {
  res = ts_transform(ap_ts, method = "none", diff = 0, seasonal_diff = 1,
                     verbose = FALSE)
  expect_equal(length(res$transformed), length(ap_ts) - frequency(ap_ts))
  expect_equal(res$params$seasonal_diff, 1L)
})

test_that("ts_transform log+diff+seasonal_diff combined", {
  res = ts_transform(ap_ts, method = "log", diff = 1, seasonal_diff = 1,
                     verbose = FALSE)
  expect_equal(res$params$method, "log")
  expect_equal(res$params$diff, 1L)
  expect_equal(res$params$seasonal_diff, 1L)
  expect_match(res$summary$steps_applied, "log")
  expect_match(res$summary$steps_applied, "diff")
  expect_match(res$summary$steps_applied, "seasonal")
})

test_that("ts_transform verbose=TRUE prints", {
  expect_output(
    ts_transform(ap_ts, method = "none", verbose = TRUE),
    "ts_transform"
  )
})

test_that("ts_transform with numeric vector and seasonal_diff auto-frequency works", {
  res = ts_transform(as.numeric(ap_ts), method = "none", seasonal_diff = 1,
                     verbose = FALSE)
  expect_equal(length(res$transformed), length(ap_ts) - 1)
  expect_equal(res$params$seasonal_diff, 1L)
})

test_that("ts_transform params stores x_vs for differencing", {
  res = ts_transform(ap_ts, method = "log", diff = 1, seasonal_diff = 0,
                     verbose = FALSE)
  expect_true(length(res$params$x_vs) > 0)
  expect_equal(as.numeric(res$params$x_vs), as.numeric(log(ap_ts)),
               tolerance = 1e-10)
})


# =============================================================================
# ts_back_transform() — Invert ts_transform()
# =============================================================================

test_that("ts_back_transform() log-only round-trip", {
  tr = ts_transform(ap_ts, method = "log", diff = 0, seasonal_diff = 0,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed)
  fc = ts_forecast(fit, h = 12)
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_s3_class(fc_orig, "tbl_df")
  expect_equal(nrow(fc_orig), 12)
  expect_true(all(c("step", "forecast") %in% names(fc_orig)))
  # forecasts should be positive (on original passenger scale)
  expect_true(all(fc_orig$forecast > 0))
})

test_that("ts_back_transform() log+diff=1 round-trip", {
  tr = ts_transform(ap_ts, method = "log", diff = 1, seasonal_diff = 0,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed)
  fc = ts_forecast(fit, h = 12)
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_s3_class(fc_orig, "tbl_df")
  expect_equal(nrow(fc_orig), 12)
  # Manual check: back-transform should roughly match log-only scale
  tr_log = ts_transform(ap_ts, method = "log", diff = 0, seasonal_diff = 0,
                        verbose = FALSE)
  fit_log = ts_sarima(tr_log$transformed)
  fc_log = ts_forecast(fit_log, h = 12)
  fc_log_orig = ts_back_transform(fc_log, tr_log$params, ap_ts)

  # diff=1 should be in the same ballpark as log-only (within 20%)
  expect_true(
    mean(abs(fc_orig$forecast - fc_log_orig$forecast) /
         fc_log_orig$forecast) < 0.20
  )
})

test_that("ts_back_transform() log+diff=1+seasonal_diff=1 round-trip", {
  tr = ts_transform(ap_ts, method = "log", diff = 1, seasonal_diff = 1,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed,
                   seasonal = list(order = c(1, 0, 1), period = 12))
  fc = ts_forecast(fit, h = 12)
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_s3_class(fc_orig, "tbl_df")
  expect_true(all(fc_orig$forecast > 0))
})

test_that("ts_back_transform() boxcox+diff=1 round-trip", {
  skip_if(!requireNamespace("forecast", quietly = TRUE),
          "forecast required")
  tr = ts_transform(ap_ts, method = "boxcox", diff = 1, seasonal_diff = 0,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed)
  fc = ts_forecast(fit, h = 12)
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_s3_class(fc_orig, "tbl_df")
  expect_true(all(fc_orig$forecast > 0))
})

test_that("ts_back_transform() preserves confidence intervals", {
  tr = ts_transform(ap_ts, method = "log", diff = 1, seasonal_diff = 0,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed)
  fc = ts_forecast(fit, h = 12, level = c(80, 95))
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_true("lo_80" %in% names(fc_orig))
  expect_true("hi_80" %in% names(fc_orig))
  expect_true("lo_95" %in% names(fc_orig))
  expect_true("hi_95" %in% names(fc_orig))
  # intervals should be monotonic: lo_95 < lo_80 < forecast < hi_80 < hi_95
  for (i in seq_len(nrow(fc_orig))) {
    expect_true(fc_orig$lo_95[i] < fc_orig$lo_80[i])
    expect_true(fc_orig$lo_80[i] < fc_orig$forecast[i])
    expect_true(fc_orig$forecast[i] < fc_orig$hi_80[i])
    expect_true(fc_orig$hi_80[i] < fc_orig$hi_95[i])
  }
})

test_that("ts_back_transform() no-diff method='none' returns identity", {
  tr = ts_transform(ap_ts, method = "none", diff = 0, seasonal_diff = 0,
                    verbose = FALSE)
  fit = ts_sarima(tr$transformed)
  fc = ts_forecast(fit, h = 5)
  fc_orig = ts_back_transform(fc, tr$params, ap_ts)

  expect_equal(fc_orig$forecast, fc$forecast, tolerance = 1e-10)
})

# =============================================================================
# Edge / error tests
# =============================================================================

test_that("ts_test errors on invalid input", {
  expect_error(ts_test("not_numeric"), "must be a numeric vector")
})

test_that("ts_ets errors on empty data frame", {
  expect_error(.validate_ts_input(data.frame()),
               "at least one numeric column")
})

test_that("ts_sarima keeps consistent lengths", {
  res = ts_sarima(ap_log)
  expect_equal(nrow(res$fitted), length(ap_log))
  expect_equal(length(res$fitted$residual), length(ap_log))
})
