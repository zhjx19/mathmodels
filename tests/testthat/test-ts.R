# =============================================================================
# Unit tests for the lightweight ts_df time-series toolkit
# =============================================================================

data(AirPassengers)
ap_ts = AirPassengers
ap_log = log(ap_ts)

# =============================================================================
# ts_df data structure
# =============================================================================

test_that("ts_df() creates a sorted two-column ts_df with metadata", {
  raw = data.frame(
    date = as.Date("2020-01-01") + c(2, 0, 1),
    sales = c(12, 10, 11)
  )

  x = ts_df(raw, time = "date", value = "sales", frequency = 7)

  expect_s3_class(x, "ts_df")
  expect_named(x, c("time", "value"))
  expect_equal(x$time, as.Date("2020-01-01") + 0:2)
  expect_equal(x$value, c(10, 11, 12))
  expect_equal(attr(x, "frequency"), 7)
  expect_equal(attr(x, "start"), as.Date("2020-01-01"))
})

test_that("ts_df() rejects duplicate times and endpoint NA values", {
  dup = data.frame(time = c(1, 1), value = c(10, 11))
  expect_error(ts_df(dup), "time must not contain duplicate values")

  leading_na = data.frame(time = 1:3, value = c(NA, 2, 3))
  expect_error(ts_df(leading_na), "First and last values must not be NA")

  trailing_na = data.frame(time = 1:3, value = c(1, 2, NA))
  expect_error(ts_df(trailing_na), "First and last values must not be NA")
})

test_that("as_ts_df() converts ts, numeric vector, and explicit data.frame inputs", {
  from_ts = as_ts_df(ts(1:4, start = c(2020, 1), frequency = 4))
  expect_s3_class(from_ts, "ts_df")
  expect_equal(from_ts$time, 1:4)
  expect_equal(from_ts$value, 1:4)
  expect_equal(attr(from_ts, "frequency"), 4)
  expect_equal(attr(from_ts, "start"), c(2020, 1))

  from_vec = as_ts_df(c(3, 5, 8), frequency = 2)
  expect_equal(from_vec$time, 1:3)
  expect_equal(from_vec$value, c(3, 5, 8))
  expect_equal(attr(from_vec, "frequency"), 2)
  expect_equal(attr(from_vec, "start"), 1)

  raw = data.frame(month = 1:3, demand = c(20, 25, 30))
  from_df = as_ts_df(raw, time = "month", value = "demand", frequency = 12, start = c(2020, 1))
  expect_equal(from_df$time, 1:3)
  expect_equal(from_df$value, c(20, 25, 30))
  expect_equal(attr(from_df, "frequency"), 12)
  expect_equal(attr(from_df, "start"), c(2020, 1))
})

test_that("as_ts_df() requires explicit columns for data.frame inputs", {
  raw = data.frame(month = 1:3, demand = c(20, 25, 30))
  expect_error(as_ts_df(raw), "time and value must be supplied")
})

test_that("is_ts_df() checks the full lightweight structure", {
  x = ts_df(data.frame(time = 1:3, value = c(10, 11, 12)))
  expect_true(is_ts_df(x))

  no_metadata = structure(
    data.frame(time = 1:3, value = c(10, 11, 12)),
    class = c("ts_df", "data.frame")
  )
  expect_false(is_ts_df(no_metadata))
})

# =============================================================================
# Completion and interpolation
# =============================================================================

test_that("complete_ts_df() fills missing time points without imputing values", {
  x = ts_df(data.frame(time = c(1, 2, 4), value = c(10, 12, 18)))

  res = complete_ts_df(x, by = 1)

  expect_s3_class(res, "ts_df")
  expect_equal(res$time, 1:4)
  expect_equal(res$value, c(10, 12, NA, 18))
  expect_equal(attr(res, "frequency"), attr(x, "frequency"))
  expect_equal(attr(res, "start"), attr(x, "start"))
})

test_that("impute_ts_df() fills internal gaps by interpolation", {
  check_case = function(values, method, expected = NULL) {
    x = ts_df(data.frame(time = c(1, 3, 5), value = values))
    res = impute_ts_df(x, by = 1, method = method)

    expect_equal(res$time, 1:5)
    expect_false(anyNA(res$value))
    expect_equal(res$value[c(1, 3, 5)], values)
    if (!is.null(expected)) {
      expect_equal(res$value, expected)
    }
  }

  check_case(c(10, 20, 30), "linear", c(10, 15, 20, 25, 30))
  check_case(c(1, 9, 25), "spline")
})

test_that("drop_na_ts_df() removes missing values and preserves metadata", {
  x = ts_df(data.frame(time = 1:4, value = c(10, NA, 12, 13)), frequency = 4)

  res = drop_na_ts_df(x)

  expect_s3_class(res, "ts_df")
  expect_equal(res$time, c(1, 3, 4))
  expect_equal(res$value, c(10, 12, 13))
  expect_equal(attr(res, "frequency"), 4)
  expect_equal(attr(res, "start"), attr(x, "start"))
})

test_that("drop_na_ts_df() removes missing values and preserves metadata", {
  x = ts_df(data.frame(time = 1:4, value = c(10, NA, 12, 13)), frequency = 4)

  res = drop_na_ts_df(x)

  expect_s3_class(res, "ts_df")
  expect_equal(res$time, c(1, 3, 4))
  expect_equal(res$value, c(10, 12, 13))
  expect_equal(attr(res, "frequency"), 4)
  expect_equal(attr(res, "start"), attr(x, "start"))
})

# =============================================================================
# Diagnostics, transformation, modeling, forecasting, and plotting
# =============================================================================

test_that("ts_stl() decomposes a seasonal ts_df into tidy components", {
  x = as_ts_df(ap_ts)

  res = ts_stl(x)

  expect_type(res, "list")
  expect_named(res, c("components", "strength", "model"))
  expect_s3_class(res$components, "tbl_df")
  expect_s3_class(res$strength, "tbl_df")
  expect_equal(nrow(res$components), nrow(x))
  expect_named(res$components, c("time", "observed", "trend", "seasonal", "remainder"))
  expect_named(res$strength, c("seasonal_strength", "trend_strength"))
})

test_that("ts_test() runs ADF and KPSS diagnostics on a clean ts_df", {
  x = as_ts_df(ap_log)

  res = ts_test(x)

  expect_s3_class(res, "tbl_df")
  expect_equal(res$test, c("ADF", "KPSS"))
  expect_type(res$statistic, "double")
  expect_type(res$p_value, "double")
})

test_that("ts_transform() returns transformed ts_df and explicit back-transform params", {
  x = as_ts_df(ap_ts)

  res = ts_transform(x, method = "log", diff = 1)

  expect_type(res, "list")
  expect_named(res, c("transformed", "params", "summary"))
  expect_s3_class(res$transformed, "ts_df")
  expect_equal(nrow(res$transformed), length(ap_ts) - 1)
  expect_equal(res$params$method, "log")
  expect_equal(res$params$diff, 1)
})

test_that("ts_back_transform() restores log-scale forecasts", {
  x = as_ts_df(ap_ts)
  res = ts_transform(x, method = "log")
  forecast_tbl = tibble::tibble(step = 1:2, forecast = log(c(500, 600)))

  restored = ts_back_transform(forecast_tbl, res$params)

  expect_equal(restored$forecast, c(500, 600), tolerance = 1e-8)
})

test_that("ts_ets() and ts_sarima() return fitted values tied to ts_df input", {
  x = as_ts_df(ap_log)

  ets = ts_ets(x)
  sarima = ts_sarima(x, stepwise = TRUE, approximation = TRUE)

  expect_setequal(names(ets), c("model_info", "parameters", "fitted", "model", "input"))
  expect_setequal(names(sarima), c("model_info", "coefficients", "fitted", "diagnostics", "model", "input"))
  expect_named(ets$fitted, c("time", "observed", "fitted", "residual"))
  expect_named(sarima$fitted, c("time", "observed", "fitted", "residual"))
  expect_s3_class(ets$input, "ts_df")
  expect_s3_class(sarima$input, "ts_df")
})

test_that("ts_arimax() fits ARIMAX with auto-selected order and xreg", {
  x = as_ts_df(ap_log)
  xreg = data.frame(trend = seq_len(nrow(x)))

  fit = ts_arimax(x, xreg = xreg, stepwise = TRUE, approximation = TRUE)

  expect_setequal(names(fit), c("model_info", "coefficients", "fitted",
                                 "diagnostics", "model", "input", "xreg"))
  expect_s3_class(fit$model, "Arima")
  expect_s3_class(fit$input, "ts_df")
  expect_true(is.matrix(fit$xreg))
  expect_equal(nrow(fit$xreg), nrow(x))
  expect_equal(fit$fitted$time, x$time)
  expect_equal(nrow(fit$coefficients), length(stats::coef(fit$model)))
})

test_that("ts_arimax() fits ARIMAX with manual order and xreg", {
  x = as_ts_df(ap_log)
  xreg = data.frame(trend = seq_len(nrow(x)))

  fit = ts_arimax(x, xreg = xreg, order = c(1, 1, 1))

  expect_s3_class(fit$model, "Arima")
  expect_equal(nrow(fit$coefficients), length(stats::coef(fit$model)))
  expect_true("xreg" %in% names(fit))
})

test_that("ts_arimax() rejects missing or mismatched xreg", {
  x = as_ts_df(ap_log)

  expect_error(ts_arimax(x), "xreg must be supplied")
  expect_error(ts_arimax(x, xreg = data.frame(a = 1:5)),
               "xreg must have .* rows")
  expect_error(ts_arimax(x, xreg = "not_a_matrix"), "xreg must be a numeric")
})

test_that("ts_forecast() works with ARIMAX model and newxreg", {
  x = as_ts_df(ap_log)
  xreg = data.frame(trend = seq_len(nrow(x)))
  fit = ts_arimax(x, xreg = xreg, order = c(1, 1, 1))

  newxreg = data.frame(trend = nrow(x) + 1:3)
  fc = ts_forecast(fit, h = 3, newxreg = newxreg)

  expect_equal(nrow(fc), 3)
  expect_named(fc, c("step", "forecast", "lo_80", "hi_80", "lo_95", "hi_95"))
})

test_that("ts_forecast() errors when newxreg is missing for ARIMAX", {
  x = as_ts_df(ap_log)
  xreg = data.frame(trend = seq_len(nrow(x)))
  fit = ts_arimax(x, xreg = xreg, order = c(1, 1, 1))

  expect_error(ts_forecast(fit, h = 3), "newxreg must be supplied")
})

test_that("ts_forecast() forecasts ETS and SARIMA results", {
  x = as_ts_df(ap_log)
  ets = ts_ets(x)
  sarima = ts_sarima(x, stepwise = TRUE, approximation = TRUE)

  ets_fc = ts_forecast(ets, h = 3, level = c(80, 95))
  sarima_fc = ts_forecast(sarima, h = 3, level = c(80, 95))

  expect_equal(nrow(ets_fc), 3)
  expect_equal(nrow(sarima_fc), 3)
  expect_named(ets_fc, c("step", "forecast", "lo_80", "hi_80", "lo_95", "hi_95"))
  expect_named(sarima_fc, c("step", "forecast", "lo_80", "hi_80", "lo_95", "hi_95"))
})

test_that("plot_ts() and plot_ts_forecast() return buildable ggplot objects", {
  x = as_ts_df(ap_log)
  fit = ts_ets(x)
  fc = ts_forecast(fit, h = 3)

  p1 = plot_ts(x)
  p2 = plot_ts_forecast(x, fc, by = 1)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p1))
  expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("plot_ts_acf() and plot_ts_pacf() return ggplot objects", {
  x = as_ts_df(ap_log)

  expect_s3_class(plot_ts_acf(x), "ggplot")
  expect_s3_class(plot_ts_pacf(x), "ggplot")
})

test_that("plot_ts_stl() returns a decomposition plot", {
  x = as_ts_df(ap_ts)
  stl_res = ts_stl(x)

  expect_s3_class(plot_ts_stl(stl_res), "patchwork")
})

test_that("plot_ts_residuals() returns a residual diagnostic plot", {
  x = as_ts_df(ap_log)
  fit = ts_sarima(x, stepwise = TRUE, approximation = TRUE)

  expect_s3_class(plot_ts_residuals(fit), "patchwork")
})
