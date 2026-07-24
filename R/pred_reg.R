# =============================================================================
# Regression prediction models for mathmodels
# Implements: LM (with stepwise), Logistic, Poisson, Negative Binomial
# Includes: model fitting, diagnostics, prediction, and visualization
# =============================================================================

# -----------------------------------------------------------------------------
# Internal Helpers
# -----------------------------------------------------------------------------

#' @importFrom stats AIC BIC as.formula binomial coef density fitted glm pchisq plogis poisson ppoints qnorm quantile residuals rnorm shapiro.test terms update vcov qt
#' @importFrom lmtest bptest dwtest
#' @importFrom car vif
#' @importFrom utils globalVariables

# Suppress R CMD check NOTE for computed ggplot2 aesthetics
globalVariables(c(".data", "x"))

# Extract variable names from a formula
.reg_formula = function(formula) {
  # response is always the LHS of the formula
  resp = as.character(attr(terms(formula), "variables"))[[2]]
  pred = attr(terms(formula), "term.labels")
  list(response = resp, predictors = pred)
}

# Validate input data
.reg_validate_data = function(data) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  if (nrow(data) < 5) stop("data must have at least 5 rows.", call. = FALSE)
}

# Compute 95% CI using model-specific variance-covariance
._ci95 = function(model, level = 0.95, se = NULL, coef_vals = NULL) {
  if (is.null(se) || is.null(coef_vals)) {
    v = vcov(model)
    se = sqrt(diag(v))
    coef_vals = coef(model)
  }
  df = model$df.residual
  z = qt((1 - level) / 2, df, lower.tail = FALSE)
  list(low = coef_vals - z * se, high = coef_vals + z * se)
}

# Build tidy coefficient tibble
._tidy_coefs = function(model, level = 0.95) {
  coefs = coef(model)
  ci = ._ci95(model, level = level)
  s = summary(model)$coefficients
  # s columns: Estimate, Std. Error, t value, Pr(>|t|)
  tibble::tibble(
    term = names(coefs),
    estimate = as.numeric(coefs),
    std.error = s[, 2],
    statistic = s[, 3],
    p.value   = s[, 4],
    conf.low  = ci$low,
    conf.high = ci$high
  )
}

# Build tidy diagnostics tibble
._tidy_diag = function(metrics) {
  tibble::tibble(
    metric = names(metrics),
    value  = round(unname(metrics), 4)
  )
}

# Generic predictor: extract response column name from fitted object
._resp_col = function(fit) {
  nm = attr(fit$model, "names")[1]
  if (nm == ".fitted") nm = names(fit$model)[1]
  nm
}

# Stepwise selection helper.
# family: NULL for lm, "binomial" for logistic, "poisson" for Poisson.
.reg_step = function(formula, data, family = NULL, direction, ...) {
  fam = if (is.null(family)) NULL else match.fun(family)()
  if (direction == "backward") {
    if (is.null(fam)) {
      model = lm(formula, data = data, ...)
    } else {
      model = glm(formula, data = data, family = fam, ...)
    }
  } else {
    f_null = update(formula, . ~ 1)
    if (is.null(fam)) {
      model = lm(f_null, data = data, ...)
    } else {
      model = glm(f_null, data = data, family = fam, ...)
    }
  }
  # Ensure data is visible inside step()'s update() calls
  model$call$data = data
  if (!is.null(fam)) model$call$family = fam
  scope = list(lower = as.formula(". ~ 1"), upper = formula)
  step(model, scope = scope, direction = direction, trace = 0)
}

# Generate new predictor data for prediction (draw from empirical distribution)
.make_newdata = function(model, model_result, n_new) {
  resp_name = attr(terms(model), "variables")[[2]]
  pred_names = setdiff(colnames(model_result$input), as.character(resp_name))
  pred_cols = model_result$input[, pred_names, drop = FALSE]
  mean_x = colMeans(pred_cols)
  sd_x = apply(pred_cols, 2, sd)
  newx = lapply(seq_along(sd_x), function(i)
    rnorm(n_new, mean = mean_x[i], sd = max(sd_x[i], 1e-8)))
  newdf = data.frame(newx)
  names(newdf) = names(mean_x)
  newdf
}


# -----------------------------------------------------------------------------
# OLS — Multivariable Linear Regression with Stepwise
# -----------------------------------------------------------------------------

#' Multivariable Linear Regression
#'
#' Fits an OLS model via `lm()`, optionally performing stepwise selection
#' using Akaike's information criterion.
#'
#' @param formula A formula, e.g. `y ~ x1 + x2`.
#' @param data A data frame containing the variables in `formula`.
#' @param step Logical. Whether to perform stepwise selection (`TRUE`/`FALSE`).
#'   Default is `FALSE`.
#' @param direction Character: `"both"` (default), `"forward"`, or `"backward"`.
#' @param ... Additional arguments passed to `lm()` and `step()`.
#'
#' @return A named list:
#' \describe{
#'   \item{model}{Fitted \code{lm} object.}
#'   \item{coefficient}{Tibble of coefficients, standard errors, t-values, p-values, 95% CIs.}
#'   \item{model_info}{Tibble of fit statistics (\code{r.squared}, \code{adj.r.squared}, AIC, BIC).}
#'   \item{residuals}{Tibble of observed, fitted, residuals.}
#'   \item{diagnostics}{Tibble of key fit metrics (\code{Breusch-Pagan}, Durbin-Watson).}
#'   \item{formula}{Original formula.}
#'   \item{input}{The original data frame.}
#' }
#'
#' @details
#' When `step = TRUE`, the function first fits a null model
#' (response only), then uses `step()` to search for the best
#' combination of terms by AIC.
#'
#' @examples
#' set.seed(42)
#' x1 = rnorm(100); x2 = rnorm(100)
#' y = 1 + 2*x1 - 3*x2 + rnorm(100, sd = 2)
#' reg_lm(y ~ x1 + x2, data = data.frame(y, x1, x2))
#'
#' @export
reg_lm = function(formula, data,
                  step = FALSE,
                  direction = c("both", "forward", "backward"),
                  ...) {
  .reg_validate_data(data)
  direction = match.arg(direction)

  if (step) {
    model = .reg_step(formula, data, family = NULL, direction, ...)
  } else {
    model = lm(formula, data = data, ...)
  }

  s = summary(model)
  n = nrow(data)
  p = length(coef(model))

  # Breusch-Pagan test
  bp_test = tryCatch(lmtest::bptest(model), error = function(e) NULL)
  bp_stat = if (!is.null(bp_test)) unname(bp_test$p.value) else NA_real_

  # Durbin-Watson test
  dw_test = tryCatch(lmtest::dwtest(model), error = function(e) NULL)
  dw_stat = if (!is.null(dw_test)) unname(dw_test$p.value) else NA_real_

  # Coefficients with CIs
  coefs_tidy = ._tidy_coefs(model)

  # Residuals tibble
  resp_col = ._resp_col(model)
  resid_tbl = tibble::tibble(
    .observed = data[[resp_col]],
    .fitted   = as.numeric(fitted(model)),
    .residual = as.numeric(residuals(model))
  )

  list(
    model       = model,
    coefficient = coefs_tidy,
    model_info  = tibble::tibble(
      r.squared     = round(s$r.squared, 4),
      adj.r.squared = round(s$adj.r.squared, 4),
      aic           = round(AIC(model), 4),
      bic           = round(BIC(model), 4)
    ),
    residuals   = resid_tbl,
    diagnostics = ._tidy_diag(c(BP_test_pvalue = bp_stat,
                                DW_test_pvalue = dw_stat)),
    formula = formula,
    input   = data
  )
}


# -----------------------------------------------------------------------------
# Logistic regression with stepwise selection
# -----------------------------------------------------------------------------

#' Logistic Regression
#'
#' Fits a binary logistic regression model via `glm(family = binomial())`,
#' with optional stepwise selection.
#'
#' @inheritParams reg_lm
#'
#' @return A named list analogous to \code{reg_lm()}, plus:
#' \describe{
#'   \item{odds_ratio}{Tibble with odds ratios and 95% CIs.}
#'   \item{residuals}{DataFrame of observed, fitted probability, residual.}
#' }
#'
#' @details
#' Same stepwise logic as \code{\link{reg_lm}}, but uses GLM with the binomial family.
#' A Hosmer-Lemeshow goodness-of-fit test is included in diagnostics.
#'
#' @examples
#' set.seed(42)
#' x1 = rnorm(100)
#' logit = 0.5 + x1
#' prob = 1 / (1 + exp(-logit))
#' y = rbinom(100, 1, prob)
#' reg_logistic(y ~ x1, data = data.frame(y, x1))
#'
#' @export
reg_logistic = function(formula, data,
                        step = FALSE,
                        direction = c("both", "forward", "backward"),
                        ...) {
  .reg_validate_data(data)
  direction = match.arg(direction)

  if (step) {
    model = .reg_step(formula, data, family = "binomial", direction, ...)
  } else {
    model = glm(formula, data = data, family = binomial(), ...)
  }

  coefs_tidy = ._tidy_coefs(model)
  # Odds ratio
  or_vals = exp(coefs_tidy$estimate)
  ci_95 = ._ci95(model)
  or_tbl = tibble::tibble(
    term       = coefs_tidy$term,
    odds_ratio = round(or_vals, 4),
    conf.low   = round(exp(ci_95$low), 4),
    conf.high  = round(exp(ci_95$high), 4)
  )

  resp_col = .reg_formula(formula)$response
  resp_raw = data[[resp_col]]
  # If factor: convert to numeric (0/1)
  if (is.factor(resp_raw)) resp_raw = as.numeric(as.character(resp_raw)) - 1L
  fitted_probs = fitted(model)
  resid_df = data.frame(
    observed         = resp_raw,
    fitted_val       = fitted_probs,
    pearson_residual = residuals(model, type = "pearson")
  )

  # Hosmer-Lemeshow (robust to duplicate quantile breaks)
  g = min(10, nrow(resid_df))
  groups = tryCatch(
    cut(fitted_probs,
        breaks = quantile(fitted_probs, seq(0, 1, length.out = g), na.rm = TRUE),
        include.lowest = TRUE),
    error = function(e) cut(fitted_probs, g, include.lowest = TRUE)
  )
  hlr = lapply(levels(groups), function(glv) {
    grp = resid_df[groups == glv, ]
    list(n = nrow(grp), obs_su = sum(grp$observed == 1),
         pred_su = sum(grp$fitted_prob))
  })
  Obs = sapply(hlr, `[[`, "obs_su")
  Exp = sapply(hlr, `[[`, "pred_su")
  HLR_stat = sum((Obs - Exp)^2 / (Exp + 1e-8) +
                 (sum(resid_df$observed) - Exp)^2 /
                 (sum(fitted_probs) + 1e-8))
  HLR_df = g - 2
  HLR_pval = pchisq(HLR_stat, HLR_df, lower.tail = FALSE)

  list(
    model       = model,
    coefficient = coefs_tidy,
    odds_ratio  = or_tbl,
    residuals   = resid_df,
    diagnostics = ._tidy_diag(c(HLR_chisq = HLR_stat,
                                HLR_df = HLR_df,
                                HLR_pvalue = HLR_pval)),
    model_info  = tibble::tibble(aic = round(AIC(model), 4),
                                  bic = round(BIC(model), 4)),
    formula     = formula,
    input       = data
  )
}


# -----------------------------------------------------------------------------
# Poisson regression
# -----------------------------------------------------------------------------

#' Poisson Regression
#'
#' Fits a count-response regression model via `glm(family = poisson())`,
#' with optional stepwise selection.
#'
#' @inheritParams reg_lm
#'
#' @return A named list analogous to \code{reg_lm()}, plus:
#' \describe{
#'   \item{dispersion}{Dispersion ratio (should be near 1).}
#'   \item{residuals}{Data frame with observed, fitted, Pearson residual.}
#' }
#'
#' @details
#' Performs an overdispersion check: dispersion > 1.5 suggests NB may be preferable.
#'
#' @examples
#' set.seed(42)
#' x1 = rnorm(100)
#' mu = exp(1 + 0.3*x1)
#' y = rpois(100, mu)
#' reg_poisson(y ~ x1, data = data.frame(y, x1))
#'
#' @export
reg_poisson = function(formula, data,
                       step = FALSE,
                       direction = c("both", "forward", "backward"),
                       ...) {
  .reg_validate_data(data)
  direction = match.arg(direction)

  if (step) {
    model = .reg_step(formula, data, family = "poisson", direction, ...)
  } else {
    model = glm(formula, data = data, family = poisson(), ...)
  }

  coefs_tidy = ._tidy_coefs(model)

  # Dispersion: sum(residual^2) / df
  pearson_resid = residuals(model, type = "pearson")
  dispersion = sum(pearson_resid^2) / model$df.residual

  resp_col = .reg_formula(formula)$response
  resid_df = data.frame(
    observed         = data[[resp_col]],
    fitted_val       = fitted(model),
    pearson_residual = pearson_resid
  )

  list(
    model       = model,
    coefficient = coefs_tidy,
    dispersion  = round(dispersion, 4),
    residuals   = resid_df,
    diagnostics = ._tidy_diag(c(dispersion = dispersion,
                                aic = AIC(model), bic = BIC(model))),
    model_info  = tibble::tibble(aic = round(AIC(model), 4),
                                  bic = round(BIC(model), 4)),
    formula     = formula,
    input       = data
  )
}


# -----------------------------------------------------------------------------
# Negative Binomial Regression
# -----------------------------------------------------------------------------

#' Negative Binomial Regression
#'
#' Fits a negative binomial regression model (for over-dispersed counts)
#' via MASS::glm.nb().
#'
#' **Note:** The `glm.nb` family does not support stepwise selection.
#' Set `step = FALSE`.
#'
#' @param formula A formula, e.g. `count ~ x1 + x2`.
#' @param data A data frame.
#' @param step Logical. Always `FALSE` for negative binomial.
#' @param ... Additional arguments passed to \code{MASS::glm.nb()}.
#'
#' @return A named list analogous to \code{reg_poisson()}, plus:
#' \describe{
#'   \item{theta}{Estimated dispersion parameter for the negative binomial distribution.}
#'   \item{residuals}{DataFrame of observed, fitted, Pearson residual.}
#' }
#'
#' @details
#' For comparison with Poisson, inspect the dispersion ratio.
#' NB is appropriate when dispersion >> 1.
#'
#' @examples
#' set.seed(42)
#' x1 = rnorm(100)
#' mu = exp(1 + 0.3*x1)
#' y = rnbinom(100, size = 1, mu = mu)
#' reg_negbin(y ~ x1, data = data.frame(y, x1))
#'
#' @export
reg_negbin = function(formula, data,
                      step = FALSE,
                      ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for reg_negbin().")
  }
  .reg_validate_data(data)
  if (step) warning("Stepwise selection is not supported for negative binomial models.")

  model = MASS::glm.nb(formula, data = data, ...)
  coefs_tidy = ._tidy_coefs(model)
  theta = model$theta

  resp_col = .reg_formula(formula)$response
  resid_df = data.frame(
    observed         = data[[resp_col]],
    fitted_val       = fitted(model),
    pearson_residual = residuals(model, type = "pearson")
  )

  list(
    model       = model,
    coefficient = coefs_tidy,
    theta       = round(theta, 4),
    residuals   = resid_df,
    diagnostics = ._tidy_diag(c(aic = AIC(model), bic = BIC(model),
                                theta = theta)),
    model_info  = tibble::tibble(aic = round(AIC(model), 4),
                                  bic = round(BIC(model), 4)),
    formula     = formula,
    input       = data
  )
}


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

#' Model Diagnostics
#'
#' Runs diagnostic tests on a fitted regression model.
#'
#' @param model_result A list returned by any \code{reg_*()} function.
#'
#' @return A named list with diagnostic tibbles. For linear models: VIF,
#'   Shapiro-Wilk normality, Breusch-Pagan, Durbin-Watson, residual stats.
#'   For logistic/poisson/negative binomial: Hosmer-Lemeshow / dispersion / deviance stats.
#' @export
reg_diagnostics = function(model_result) {
  m = model_result$model
  res = list()

  if (inherits(m, "glm")) {
    # GLM diagnostics
    family_type = if (!is.null(m$family$family)) {
      tolower(as.character(m$family$family))
    } else {
      "unknown"
    }

    if (family_type == "binomial") {
      hlr_pval = tryCatch(
        model_result$diagnostics$value[
          model_result$diagnostics$metric == "HLR_pvalue"],
        error = function(e) NA_real_
      )
      res$hosmer_lemeshow = tibble::tibble(
        chi_sq = NA_real_, df = NA_integer_,
        p.value = round(hlr_pval, 6))
    } else {
      dispersion = tryCatch(
        sum(residuals(m, type = "pearson")^2) / m$df.residual,
        error = function(e) NA_real_)
      res$dispersion = tibble::tibble(value = round(dispersion, 4))
    }

    res$deviance_residual_stats = tibble::tibble(
      statistic = c("mean", "sd", "min", "max"),
      value = round(c(
        mean(residuals(m)),
        sd(residuals(m)),
        min(residuals(m)),
        max(residuals(m))
      ), 4)
    )

  } else if (inherits(m, "lm")) {
    # VIF
    car_vif = tryCatch(car::vif(m), error = function(e) NULL)
    if (!is.null(car_vif)) {
      res$vif = tibble::tibble(term = names(car_vif), vif = round(car_vif, 4))
    } else {
      res$vif = tibble::tibble(term = character(0), vif = numeric(0))
    }

    # Shapiro-Wilk
    resid_vec = as.numeric(residuals(m))
    sw = shapiro.test(resid_vec)
    res$shapiro_wilk = tibble::tibble(
      statistic = round(sw$statistic, 4),
      p.value   = round(sw$p.value, 6))

    # Residual stats
    res$residual_stats = tibble::tibble(
      statistic = c("mean", "sd", "min", "max"),
      value = round(c(mean(resid_vec), sd(resid_vec),
                      min(resid_vec), max(resid_vec)), 4)
    )

    # Breusch-Pagan
    res$breusch_pagan = tryCatch({
      tbl = lmtest::bptest(m)
      tibble::tibble(statistic = round(tbl$statistic, 4),
                     df = tbl$df, p.value = round(tbl$p.value, 6))
    }, error = function(e) tibble::tibble(
      statistic = NA_real_, df = NA_integer_, p.value = NA_real_))

    # Durbin-Watson
    res$dw_test = tryCatch({
      tbl = lmtest::dwtest(m)
      tibble::tibble(statistic = round(tbl$statistic, 4),
                     p.value = round(tbl$p.value, 6))
    }, error = function(e) tibble::tibble(
      statistic = NA_real_, p.value = NA_real_))
  }

  res
}


# -----------------------------------------------------------------------------
# Prediction
# -----------------------------------------------------------------------------

#' Predictions
#'
#' Generates point predictions and confidence/ prediction intervals for new data.
#'
#' @param model_result A fitted result from any \code{reg_*()} function.
#' @param n_new Number of new observations to predict. Random values drawn from
#'   the fitted model's residuals are used for predictor data.
#' @param level Confidence/interval level. Defaults to 0.95.
#'
#' @return A named list with:
#' \describe{
#'   \item{predictions}{Tibble with predicted values.}
#'   \item{confidence_interval}{Tibble with lower and upper bounds.}
#' }
#' @export
reg_predict = function(model_result, n_new = 10, level = 0.95) {
  m = model_result$model

  if (inherits(m, "glm")) {
    fam = m$family
    if (fam$link == "logit") {
      return(.predict_logistic(m, model_result, n_new = n_new, level = level))
    } else if (fam$link == "log") {
      return(.predict_count(m, model_result, n_new = n_new, level = level))
    }
    stop("Unrecognized GLM family.", call. = FALSE)
  } else if (inherits(m, "lm")) {
    return(.predict_lm(m, model_result, n_new = n_new, level = level))
  } else {
    stop("Unsupported model type.", call. = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Prediction helpers
# -----------------------------------------------------------------------------

#' @importFrom purrr map map_dbl map_dfr
.predict_lm = function(model, model_result, n_new, level) {
  newdf = .make_newdata(model, model_result, n_new)
  pred = predict(model, newdata = newdf, interval = "prediction", level = level)

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_new),
      .predict = pred[, "fit"]),
    confidence_interval = tibble::tibble(
      lower = round(as.numeric(pred[, "lwr"]), 4),
      upper = round(as.numeric(pred[, "upr"]), 4))
  )
}

#' @importFrom purrr map map_dbl map_dfr
.predict_logistic = function(model, model_result, n_new, level) {
  newdf = .make_newdata(model, model_result, n_new)
  pred = predict(model, newdata = newdf, type = "response", se.fit = TRUE)
  z = qnorm(1 - (1 - level) / 2)
  ci_lower = plogis(log(pred$fit / (1 - pred$fit + 1e-8)) - z * pred$se.fit)
  ci_upper = plogis(log(pred$fit / (1 - pred$fit + 1e-8)) + z * pred$se.fit)

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_new),
      .predict = round(pred$fit, 4)),
    confidence_interval = tibble::tibble(
      lower = round(ci_lower, 4),
      upper = round(ci_upper, 4))
  )
}

#' @importFrom purrr map map_dbl map_dfr
.predict_count = function(model, model_result, n_new, level) {
  newdf = .make_newdata(model, model_result, n_new)
  pred = predict(model, newdata = newdf, type = "response", se.fit = TRUE)
  z = qnorm(1 - (1 - level) / 2)
  ci_lower = pmax(pred$fit - z * pred$se.fit, 0)
  ci_upper = pred$fit + z * pred$se.fit

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_new),
      .predict = round(pred$fit, 4)),
    confidence_interval = tibble::tibble(
      lower = round(as.numeric(ci_lower), 4),
      upper = round(as.numeric(ci_upper), 4))
  )
}


# -----------------------------------------------------------------------------
# Visualization helpers
# -----------------------------------------------------------------------------

.resid_plot_lm = function(model_result) {
  resid = model_result$residuals
  n_rows = nrow(resid)

  p1 = ggplot2::ggplot(resid, ggplot2::aes(x = .data$.fitted, y = .data$.residual)) +
    ggplot2::geom_point(size = 2, alpha = 0.6, colour = "#2C3E50") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#BDC3C7") +
    ggplot2::labs(title = "Residuals vs Fitted",
                  x = "Fitted Values", y = "Residuals")

  # Q-Q
  qq_data = tibble::tibble(
    theoretical = stats::qnorm(ppoints(n_rows)),
    sample = sort(resid$.residual)
  )
  p2 = ggplot2::ggplot(qq_data, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_point(colour = "#2C3E50", size = 1.5, alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#E74C3C") +
    ggplot2::labs(title = "Q-Q Plot", x = "Theoretical", y = "Sample Quantiles")

  # Histogram
  p3 = ggplot2::ggplot(resid, ggplot2::aes(x = .data$.residual)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 30, fill = "#2C3E50", colour = "white", alpha = 0.7) +
    ggplot2::geom_density(colour = "#E74C3C") +
    ggplot2::labs(title = "Residual Histogram", x = "Residuals", y = "Density")

  # ACF
  resid_clean = resid$.residual[!is.na(resid$.residual)]
  acf_obj = if (length(resid_clean) > 2) {
    stats::acf(resid_clean, plot = FALSE)
  } else {
    NULL
  }
  ci = stats::qnorm(0.975) / sqrt(n_rows)
  acf_data = tibble::tibble(
    lag = as.numeric(acf_obj$lag),
    acf = as.numeric(acf_obj$acf),
    ci  = ci
  )
  p4 = ggplot2::ggplot(acf_data, ggplot2::aes(x = .data$lag, y = .data$acf)) +
    ggplot2::geom_hline(yintercept = c(-ci, ci), linetype = "dashed", colour = "#95A5A6") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0),
                          colour = "#2C3E50", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#2C3E50") +
    ggplot2::labs(title = "ACF of Residuals", x = "Lag", y = "ACF")

  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2) +
    patchwork::plot_annotation(
      title = "Linear Model Diagnostics",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
}

.resid_plot_glm = function(model_result) {
  resid = model_result$residuals
  n_rows = nrow(resid)

  p1 = ggplot2::ggplot(resid, ggplot2::aes(x = .data$fitted_val, y = .data$pearson_residual)) +
    ggplot2::geom_point(size = 2, alpha = 0.6, colour = "#2C3E50") +
    ggplot2::geom_hline(yintercept = c(-2, 0, 2), linetype = "dashed", colour = "#BDC3C7") +
    ggplot2::labs(title = "Pearson Residuals vs Fitted",
                  x = "Fitted", y = "Pearson Resid.")

  # QQ for GLM
  norm_resid = sort(stats::rstandard(model_result$model))
  n_norm = length(norm_resid)
  qq_data = tibble::tibble(
    theoretical = stats::qnorm(ppoints(n_norm)),
    sample = norm_resid
  )
  p2 = ggplot2::ggplot(qq_data, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_point(colour = "#2C3E50", size = 1.5, alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#E74C3C") +
    ggplot2::labs(title = "QQ Plot", x = "Theoretical", y = "Sample Quantiles")

  patchwork::wrap_plots(p1, p2, ncol = 2) +
    patchwork::plot_annotation(
      title = paste(model_result$coefficients$term[1], "Model Diagnostics"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
}

#' Residual Diagnostics
#'
#' Returns a multi-panel diagnostic plot from a fitted regression model.
#'
#' @param model_result A list returned by any \code{reg_*()} function.
#'
#' @return A ggplot-based composite.
#' @export
reg_plot_residuals = function(model_result) {
  if (inherits(model_result$model, "glm")) {
    return(.resid_plot_glm(model_result))
  } else if (inherits(model_result$model, "lm")) {
    return(.resid_plot_lm(model_result))
  } else {
    stop("Unsupported model for residual plot.")
  }
}


# -----------------------------------------------------------------------------
# Prediction plots
# -----------------------------------------------------------------------------

#' @importFrom stats predict
#' @keywords internal
.prediction_plot_generic = function(model_result, type = c("fitted", "predicted"),
                                    show_ci = TRUE) {
  type = match.arg(type)
  m = model_result$model

  if (inherits(m, "glm")) {
    fam = m$family
    form = if (is.null(m$formula)) model_result$formula else m$formula
    pred_names = setdiff(colnames(m$model), .reg_formula(form)$response)
    x_var = pred_names[1]
    x_data = m$model[[x_var]]
    fitted_vals = fitted(m)

    if (fam$link == "logit") {
      p = ggplot2::ggplot(
        data.frame(x = x_data, fitted = fitted_vals),
        ggplot2::aes(x = x, y = fitted)) +
        ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
        ggplot2::labs(title = "Predicted Probabilities", x = x_var, y = "P(y=1)")
    } else {
      p = ggplot2::ggplot(
        data.frame(x = x_data, fitted = fitted_vals),
        ggplot2::aes(x = x, y = fitted)) +
        ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
        ggplot2::labs(title = "Fitted Count Values", x = x_var, y = "Count")
    }
  } else if (inherits(m, "lm")) {
    x_var = colnames(m$model)[2]
    x_data = m$model[, x_var]
    fitted_vals = fitted(m)

    if (show_ci) {
      ci = predict(m, interval = "confidence")
      p = ggplot2::ggplot(
        data.frame(x = x_data, fitted = fitted_vals,
                   ci_lo = ci[, "lwr"], ci_hi = ci[, "upr"]),
        ggplot2::aes(x = x, y = fitted)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi),
                             alpha = 0.15, fill = "#2980B9") +
        ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
        ggplot2::labs(title = "Fitted Values with 95% CI",
                      x = x_var, y = "Response")
    } else {
      p = ggplot2::ggplot(
        data.frame(x = x_data, fitted = fitted_vals),
        ggplot2::aes(x = x, y = fitted)) +
        ggplot2::geom_line(colour = "#2C3E50", linewidth = 1) +
        ggplot2::labs(title = "Fitted Values", x = x_var, y = "Response")
    }
  } else {
    stop("Unsupported model for prediction plot.")
  }

  p
}

#' Prediction Plot
#'
#' Returns a fitted-vs-observed line plot, with optional confidence bands.
#'
#' @param model_result A fitted result from any \code{reg_*()} function.
#' @param type Character. Currently supports `"fitted"` only.
#' @param show_ci Logical. Whether to overlay confidence intervals. Default `TRUE`.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#' @export
reg_plot_predict = function(model_result,
                             type = c("fitted", "predicted"),
                             show_ci = TRUE) {
  .prediction_plot_generic(model_result, type, show_ci = show_ci)
}
