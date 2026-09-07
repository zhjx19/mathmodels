# =============================================================================
# Unit tests for the regression prediction module pred_reg.R
# =============================================================================

# --- Regression data helpers ---

test_that("reg_formula() extracts variable names from a formula", {
  f = y ~ x1 + x2 + x3

  expect_equal(mathmodels:::.reg_formula(f)$response, "y")
  expect_equal(mathmodels:::.reg_formula(f)$predictors, c("x1", "x2", "x3"))

  # Single predictor
  f2 = response ~ factor
  expect_equal(mathmodels:::.reg_formula(f2)$predictors, "factor")
})

test_that(".reg_formula() works with interaction terms", {
  f = y ~ x1 * x2

  expect_equal(mathmodels:::.reg_formula(f)$response, "y")
  expect_setequal(mathmodels:::.reg_formula(f)$predictors, c("x1", "x2", "x1:x2"))
})

test_that(".reg_validate_data() checks required arguments", {
  df = data.frame(a = 1:10, b = 1:10)

  # Should work
  expect_no_error(mathmodels:::.reg_validate_data(df))

  # Empty data.frame
  expect_error(mathmodels:::.reg_validate_data(data.frame()))
})

# --- Regression fit functions ---

test_that("reg_lm() returns a complete result list", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  res = reg_lm(y ~ x1 + x2, data = df)

  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "diagnostics", "formula", "input"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_s3_class(res$coefficient, "tbl_df")
  expect_s3_class(res$diagnostics, "tbl_df")
  expect_s3_class(res$residuals, "tbl_df")
  expect_equal(nrow(res$residuals), nrow(df))
})

test_that("reg_lm() coefficients contain expected columns", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  res = reg_lm(y ~ x1 + x2, data = df)

  expect_setequal(names(res$coefficient), c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
})

test_that("reg_lm() model_info contains fit statistics", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  res = reg_lm(y ~ x1, data = df)

  expect_setequal(names(res$model_info), c("r.squared", "adj.r.squared", "aic", "bic"))
  expect_true(is.numeric(res$model_info$r.squared[1]))
})

test_that("reg_lm() diagnostics contain BP and DW test p-values", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  res = reg_lm(y ~ x1, data = df)

  expect_true("BP_test_pvalue" %in% res$diagnostics$metric)
  expect_true("DW_test_pvalue" %in% res$diagnostics$metric)
})

# --- Logistic regression ---

test_that("reg_logistic() returns valid model results", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x,
    z = rnorm(n)
  )

  res = reg_logistic(y ~ x + z, data = df)

  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "odds_ratio", "residuals", "diagnostics", "model_info", "formula", "input"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_s3_class(res$coefficient, "tbl_df")
  expect_s3_class(res$odds_ratio, "tbl_df")
})

test_that("reg_logistic() includes odds ratio with 95% CI", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x
  )

  res = reg_logistic(y ~ x, data = df)

  expect_equal(nrow(res$odds_ratio), nrow(res$coefficient))
  expect_setequal(names(res$odds_ratio), c("term", "odds_ratio", "conf.low", "conf.high"))
})

# --- Poisson regression ---

test_that("reg_poisson() returns valid model results", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  mu = exp(1 + 0.5 * x)
  df = data.frame(
    count = rpois(n, mu),
    x = x
  )

  res = reg_poisson(count ~ x, data = df)

  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "dispersion", "residuals", "diagnostics", "model_info", "formula", "input"))
  expect_s3_class(res$coefficient, "tbl_df")
  expect_true(is.numeric(res$dispersion))
})

test_that("reg_poisson() reports dispersion statistic", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  mu = exp(1 + 0.5 * x)
  df = data.frame(
    count = rpois(n, mu),
    x = x
  )

  res = reg_poisson(count ~ x, data = df)

  expect_true(any(c("dispersion") %in% res$diagnostics$metric))
})

# --- Negative binomial regression ---

test_that("reg_negbin() returns valid model results", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  mu = exp(1 + 0.5 * x)
  df = data.frame(
    count = rnbinom(n, size = 1, mu = mu),
    x = x
  )

  res = reg_negbin(count ~ x, data = df)

  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "theta", "residuals", "diagnostics", "model_info", "formula", "input"))
  expect_s3_class(res$coefficient, "tbl_df")
  expect_true(is.numeric(res$theta))
})

# --- Prediction ---

test_that("reg_predict() produces predictions matching model type", {
  df = data.frame(
    y = rnorm(30),
    x1 = rnorm(30),
    x2 = rnorm(30)
  )

  fit = reg_lm(y ~ x1 + x2, data = df)
  preds = reg_predict(fit, n_new = 5)

  expect_type(preds, "list")
  expect_named(preds, c("predictions", "confidence_interval", "prediction_interval", "newdata"))
  expect_s3_class(preds$predictions, "tbl_df")
  expect_equal(nrow(preds$predictions), 5)
  expect_named(preds$confidence_interval, c("lower", "upper"))
  expect_equal(length(preds$confidence_interval$lower), 5)
})

test_that("reg_predict() works for logistic regression", {
  set.seed(42)
  n = 100
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x
  )

  fit = reg_logistic(y ~ x, data = df)
  preds = reg_predict(fit, n_new = 3)

  expect_type(preds, "list")
  expect_named(preds, c("predictions", "confidence_interval", "newdata"))
  expect_s3_class(preds$predictions, "tbl_df")
  expect_equal(nrow(preds$predictions), 3)
  # For logistic, probabilities should be in [0, 1]
  expect_true(all(preds$predictions$.predict >= 0 & preds$predictions$.predict <= 1))
})

# --- Residual plot edge cases ---

test_that("plot_reg_residuals() works for GLM", {
  set.seed(42)
  n = 50
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x
  )

  fit = reg_logistic(y ~ x, data = df)
  p = plot_reg_residuals(fit)
  expect_s3_class(p, "patchwork")
})

test_that("reg_predict() error handling for missing model type", {
  fake = list(model_info = tibble::tibble(model_type = "unknown"))
  expect_error(reg_predict(fake), "Unsupported model type")
})

# --- Diagnostics ---

test_that("reg_diagnostics() returns valid diagnostic summary for lm", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  fit = reg_lm(y ~ x1 + x2, data = df)
  diag = reg_diagnostics(fit)

  expect_type(diag, "list")
  expect_setequal(names(diag), c("vif", "shapiro_wilk", "breusch_pagan", "dw_test", "residual_stats"))
  expect_s3_class(diag$vif, "tbl_df")
  expect_s3_class(diag$shapiro_wilk, "tbl_df")
  expect_s3_class(diag$residual_stats, "tbl_df")
})

test_that("reg_diagnostics() vif returns tibble with v values", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  fit = reg_lm(y ~ x1 + x2, data = df)
  diag = reg_diagnostics(fit)

  expect_s3_class(diag$vif, "tbl_df")
  expect_setequal(names(diag$vif), c("term", "vif"))
})

test_that("reg_diagnostics() shapiro_wilk returns valid p-value", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  fit = reg_lm(y ~ x1, data = df)
  diag = reg_diagnostics(fit)

  sw_pval = diag$shapiro_wilk$p.value[1]
  expect_true(is.numeric(sw_pval) && sw_pval >= 0 && sw_pval <= 1)
})

test_that("reg_diagnostics() residual_stats includes basic stats", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50)
  )

  fit = reg_lm(y ~ x1, data = df)
  diag = reg_diagnostics(fit)

  rs = diag$residual_stats
  expect_setequal(rs$statistic, c("mean", "sd", "min", "max"))
  expect_equal(nrow(rs), 4)
})

test_that("reg_diagnostics() returns valid diagnostics for logistic model", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x
  )

  fit = reg_logistic(y ~ x, data = df)
  diag = reg_diagnostics(fit)

  expect_type(diag, "list")
  # GLM diagnostics should include deviance stats
  expect_true(any(c("hosmer_lemeshow", "dispersion") %in% names(diag)))
})

# --- Plotting ---

test_that("plot_reg_residuals() returns a ggplot object", {
  df = data.frame(
    y = rnorm(50),
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  fit = reg_lm(y ~ x1 + x2, data = df)

  p = plot_reg_residuals(fit)
  expect_s3_class(p, "ggplot")
})

test_that("plot_reg_residuals() works for Poisson model", {
  set.seed(42)
  n = 100
  x = rnorm(n)
  mu = exp(1 + 0.3 * x)
  df = data.frame(
    count = rpois(n, mu),
    x = x
  )

  fit = reg_poisson(count ~ x, data = df)
  p = plot_reg_residuals(fit)

  expect_s3_class(p, "patchwork")
})

test_that("plot_reg_residuals() works for GLM returns patchwork", {
  set.seed(42)
  n = 50
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(
    y = as.factor(rbinom(n, 1, prob)),
    x = x
  )

  fit = reg_logistic(y ~ x, data = df)
  p = plot_reg_residuals(fit)

  # GLM residual plot returns patchwork (2-panel layout)
  expect_s3_class(p, "patchwork")
})

# --- Value-level regression fixes ---

test_that("reg_logistic() HL statistic is sane and matches reference implementation", {
  set.seed(123)
  n = 200
  x = rnorm(n)
  prob = plogis(0.3 + 0.8 * x)
  y = rbinom(n, 1, prob)
  df = data.frame(y = y, x = x)

  res = reg_logistic(y ~ x, data = df)

  hl_chisq = res$diagnostics$value[res$diagnostics$metric == "HLR_chisq"]
  hl_df = res$diagnostics$value[res$diagnostics$metric == "HLR_df"]
  hl_p = res$diagnostics$value[res$diagnostics$metric == "HLR_pvalue"]

  expect_true(is.finite(hl_chisq))
  expect_gt(hl_chisq, 0)
  expect_lt(hl_chisq, 50)
  expect_gt(hl_p, 0)
  expect_lt(hl_p, 1)

  # Hand-rolled reference: standard HL with decile grouping of fitted probs
  # (g = 10 groups require g + 1 = 11 quantile breaks)
  p_hat = as.numeric(fitted(res$model))
  brks = quantile(p_hat, seq(0, 1, length.out = 11))
  grp = cut(p_hat, breaks = brks, include.lowest = TRUE, labels = FALSE)
  stat_ref = 0
  n_groups = 0
  for (k in sort(unique(grp))) {
    sel = grp == k
    o = sum(y[sel])
    e = sum(p_hat[sel])
    ng = sum(sel)
    if (e > 0 && (ng - e) > 0) {
      stat_ref = stat_ref + (o - e)^2 / (e * (1 - e / ng))
      n_groups = n_groups + 1
    }
  }
  expect_equal(hl_chisq, stat_ref, tolerance = 1e-3)
  expect_equal(hl_df, n_groups - 2)
})

test_that("reg_logistic() handles 'no'/'yes' factor response correctly", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  y01 = rbinom(n, 1, prob)
  df = data.frame(
    y = factor(ifelse(y01 == 1, "yes", "no"), levels = c("no", "yes")),
    x = x
  )

  res = reg_logistic(y ~ x, data = df)

  # observed must be 0/1 with 1 = event ("yes"), not NA
  expect_setequal(unique(res$residuals$observed), c(0, 1))
  expect_equal(mean(res$residuals$observed), mean(y01), tolerance = 1e-8)
  expect_true(all(res$residuals$fitted_val >= 0 & res$residuals$fitted_val <= 1))
  expect_true(is.finite(res$diagnostics$value[res$diagnostics$metric == "HLR_chisq"]))

  preds = reg_predict(res, n_new = 5)
  expect_true(all(preds$predictions$.predict >= 0 & preds$predictions$.predict <= 1))
})

test_that("reg_logistic() prediction CI has plausible width", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(y = factor(rbinom(n, 1, prob)), x = x)

  fit = reg_logistic(y ~ x, data = df)
  preds = reg_predict(fit, n_new = 10)

  lower = preds$confidence_interval$lower
  upper = preds$confidence_interval$upper
  expect_true(all(lower < preds$predictions$.predict))
  expect_true(all(preds$predictions$.predict < upper))
  # The old response-scale-SE bug gave widths around 1e-4
  expect_true(all((upper - lower) > 0.02))
  expect_true(all(lower >= 0 & upper <= 1))
})

test_that("reg_predict() lm returns both confidence and prediction intervals", {
  set.seed(42)
  df = data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50))

  fit = reg_lm(y ~ x1 + x2, data = df)
  preds = reg_predict(fit, n_new = 5)

  expect_named(preds, c("predictions", "confidence_interval", "prediction_interval", "newdata"))
  expect_equal(nrow(preds$confidence_interval), 5)
  expect_equal(nrow(preds$prediction_interval), 5)
  # Prediction interval must be wider than the confidence interval
  expect_true(all(preds$prediction_interval$lower < preds$confidence_interval$lower))
  expect_true(all(preds$prediction_interval$upper > preds$confidence_interval$upper))
})

test_that("reg_predict() uses supplied newdata and returns it", {
  set.seed(42)
  df = data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50))

  fit = reg_lm(y ~ x1 + x2, data = df)
  nd = data.frame(x1 = c(1, 2), x2 = c(-1, 0))

  preds = reg_predict(fit, newdata = nd)

  expect_equal(nrow(preds$predictions), 2)
  expect_equal(nrow(preds$confidence_interval), 2)
  expect_equal(preds$newdata$x1, nd$x1)
  expect_equal(preds$newdata$x2, nd$x2)
  manual = as.numeric(predict(fit$model, newdata = nd))
  expect_equal(as.numeric(preds$predictions$.predict), manual, tolerance = 1e-4)
})

test_that("reg_predict() handles factor predictors without error", {
  set.seed(42)
  n = 60
  g = factor(sample(c("A", "B"), n, replace = TRUE))
  df = data.frame(y = rnorm(n) + ifelse(g == "B", 1, 0), g = g)

  fit = reg_lm(y ~ g, data = df)

  expect_no_error(preds <- reg_predict(fit, n_new = 4))
  expect_equal(nrow(preds$predictions), 4)
  expect_type(preds$newdata$g, "integer")  # factor preserved
})

test_that("reg_lm() handles rank-deficient design (aliased coefficients)", {
  set.seed(42)
  a = rnorm(50)
  b = 2 * a + rnorm(50, sd = 0.1)
  c = a + b
  df = data.frame(y = rnorm(50), a = a, b = b, c = c)

  expect_no_error(res <- reg_lm(y ~ a + b + c, data = df))

  # Aliased term carries NA in coef() but must be dropped from the tidy table
  n_aliased = sum(is.na(coef(res$model)))
  expect_equal(nrow(res$coefficient), length(coef(res$model)) - n_aliased)
  expect_true(!any(is.na(res$coefficient$estimate)))
  expect_true(!any(is.na(res$coefficient$std.error)))
})

test_that("reg_diagnostics() handles n > 5000 residuals for shapiro.test", {
  set.seed(42)
  n = 5001
  df = data.frame(y = rnorm(n), x1 = rnorm(n))

  fit = reg_lm(y ~ x1, data = df)

  expect_no_error(diag <- reg_diagnostics(fit))
  expect_s3_class(diag$shapiro_wilk, "tbl_df")
})

test_that("reg_diagnostics() recomputes Hosmer-Lemeshow from stored predictions", {
  set.seed(42)
  n = 200
  x = rnorm(n)
  prob = plogis(0.5 + x)
  df = data.frame(y = as.factor(rbinom(n, 1, prob)), x = x)

  fit = reg_logistic(y ~ x, data = df)
  diag = reg_diagnostics(fit)

  hl = diag$hosmer_lemeshow
  expect_s3_class(hl, "tbl_df")
  expect_setequal(names(hl), c("chi_sq", "df", "p.value"))
  expect_true(is.finite(hl$chi_sq))
  expect_true(hl$chi_sq > 0 && hl$chi_sq < 50)
  expect_equal(hl$df, 8)
  expect_gt(hl$p.value, 0)
  expect_lt(hl$p.value, 1)
})

test_that("reg_lm() handles transformed response (log(y))", {
  set.seed(42)
  n = 50
  x = rnorm(n)
  y = exp(1 + 0.5 * x + rnorm(n, sd = 0.3))
  df = data.frame(x = x, y = y)

  res = reg_lm(log(y) ~ x, data = df)

  expect_equal(nrow(res$residuals), n)
  expect_setequal(names(res$residuals), c(".observed", ".fitted", ".residual"))
  expect_true(!any(is.na(res$residuals$.observed)))
})

test_that("residuals output is a tibble for all model types", {
  set.seed(42)
  n = 100
  x = rnorm(n)

  lm_res = reg_lm(y ~ x, data = data.frame(y = rnorm(n), x = x))
  logistic_res = reg_logistic(
    y ~ x, data = data.frame(y = as.factor(rbinom(n, 1, plogis(x))), x = x))
  poisson_res = reg_poisson(
    count ~ x, data = data.frame(count = rpois(n, exp(1 + 0.3 * x)), x = x))
  negbin_res = reg_negbin(
    count ~ x, data = data.frame(count = rnbinom(n, size = 1, mu = exp(1 + 0.3 * x)), x = x))

  expect_s3_class(lm_res$residuals, "tbl_df")
  expect_s3_class(logistic_res$residuals, "tbl_df")
  expect_s3_class(poisson_res$residuals, "tbl_df")
  expect_s3_class(negbin_res$residuals, "tbl_df")
})
