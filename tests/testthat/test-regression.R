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
  expect_named(preds, c("predictions", "confidence_interval"))
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
  expect_named(preds, c("predictions", "confidence_interval"))
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
