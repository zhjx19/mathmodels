# =============================================================================
# Regression prediction models for mathmodels
# Implements: LM (with stepwise), Logistic, Poisson, Negative Binomial
# Includes: model fitting, diagnostics, prediction, and visualization
# =============================================================================

# -----------------------------------------------------------------------------
# Internal Helpers
# -----------------------------------------------------------------------------

#' @importFrom stats AIC BIC as.formula binomial coef density fitted glm pchisq plogis poisson ppoints qnorm quantile residuals rnorm shapiro.test step terms update vcov qt
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

# Build tidy coefficient tibble.
# Built from summary()$coefficients so aliased (rank-deficient) terms,
# which only carry NA in coef(), are excluded consistently.
._tidy_coefs = function(model, level = 0.95) {
  s = summary(model)$coefficients
  # s columns: Estimate, Std. Error, t/z value, Pr(>|t/z|)
  est = as.numeric(s[, 1])
  se  = as.numeric(s[, 2])
  crit = qt((1 - level) / 2, model$df.residual, lower.tail = FALSE)
  tibble::tibble(
    term      = rownames(s),
    estimate  = est,
    std.error = se,
    statistic = as.numeric(s[, 3]),
    p.value   = as.numeric(s[, 4]),
    conf.low  = est - crit * se,
    conf.high = est + crit * se
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

# Generate new predictor data for prediction.
# When newdata is supplied it is used directly; otherwise rows are drawn from
# the empirical distributions of the predictors: normal approximation
# (mean/sd) for numeric columns, levels sampled with observed proportions
# for factor/character columns.
.make_newdata = function(model, model_result, n_new, newdata = NULL) {
  if (!is.null(newdata)) return(as.data.frame(newdata))
  resp_name = attr(terms(model), "variables")[[2]]
  pred_names = setdiff(colnames(model_result$input), as.character(resp_name))
  pred_cols = model_result$input[, pred_names, drop = FALSE]
  newdf = as.data.frame(lapply(pred_cols, function(col) {
    if (is.factor(col) || is.character(col)) {
      tab = table(col)
      out = sample(names(tab), size = n_new, replace = TRUE,
                   prob = as.numeric(tab))
      if (is.factor(col)) factor(out, levels = levels(col)) else out
    } else {
      rnorm(n_new, mean = mean(col, na.rm = TRUE),
            sd = max(sd(col, na.rm = TRUE), 1e-8))
    }
  }), stringsAsFactors = FALSE)
  newdf
}

# Hosmer-Lemeshow goodness-of-fit test for binary logistic models.
# observed: 0/1 vector (1 = event); fitted_probs: predicted probabilities.
# Groups observations by deciles of predicted probability (g = 10) and
# computes the standard statistic sum_g (O_g - E_g)^2 / (E_g (1 - E_g/n_g)).
.hosmer_lemeshow = function(observed, fitted_probs, g = 10) {
  observed = as.numeric(observed)
  fitted_probs = as.numeric(fitted_probs)
  if (length(fitted_probs) < 5 || anyNA(fitted_probs) || anyNA(observed)) {
    return(list(stat = NA_real_, df = NA_integer_, p.value = NA_real_,
                n_groups = 0L))
  }
  g = min(g, length(fitted_probs))
  grp_idx = tryCatch(
    cut(fitted_probs,
        breaks = quantile(fitted_probs, seq(0, 1, length.out = g + 1), na.rm = TRUE),
        include.lowest = TRUE, labels = FALSE),
    error = function(e) cut(fitted_probs, g, include.lowest = TRUE, labels = FALSE)
  )
  hl_groups = lapply(seq_len(g), function(k) {
    sel = which(grp_idx == k)
    if (length(sel) == 0) return(NULL)
    list(n = length(sel),
         obs = sum(observed[sel] == 1),
         exp = sum(fitted_probs[sel]))
  })
  hl_groups = Filter(Negate(is.null), hl_groups)
  # Skip groups with a degenerate denominator (E_g = 0 or n_g - E_g = 0)
  hl_groups = Filter(function(x) x$exp > 0 && (x$n - x$exp) > 0, hl_groups)
  n_groups = length(hl_groups)
  if (n_groups < 3) {
    return(list(stat = NA_real_, df = NA_integer_, p.value = NA_real_,
                n_groups = n_groups))
  }
  stat = sum(vapply(hl_groups, function(x)
    (x$obs - x$exp)^2 / (x$exp * (1 - x$exp / x$n)), numeric(1)))
  df = n_groups - 2L
  list(stat = stat, df = df,
       p.value = pchisq(stat, df, lower.tail = FALSE), n_groups = n_groups)
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
#' res = reg_lm(y ~ x1 + x2, data = data.frame(y, x1, x2))
#' res$coefficient
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
  obs_vec = data[[resp_col]]
  if (is.null(obs_vec)) obs_vec = as.numeric(stats::model.response(model$model))
  resid_tbl = tibble::tibble(
    .observed = obs_vec,
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
#'   \item{residuals}{Tibble of observed (0/1), fitted probability, Pearson residual.}
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
#' res = reg_logistic(y ~ x1, data = data.frame(y, x1))
#' res$coefficient
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
  # Map the response to 0/1 (glm convention: success = second factor level)
  if (is.factor(resp_raw)) resp_raw = as.integer(resp_raw == levels(resp_raw)[2])
  fitted_probs = as.numeric(fitted(model))
  resid_df = tibble::tibble(
    observed         = resp_raw,
    fitted_val       = fitted_probs,
    pearson_residual = as.numeric(residuals(model, type = "pearson"))
  )

  # Hosmer-Lemeshow goodness-of-fit test
  hl = .hosmer_lemeshow(resid_df$observed, fitted_probs)

  list(
    model       = model,
    coefficient = coefs_tidy,
    odds_ratio  = or_tbl,
    residuals   = resid_df,
    diagnostics = ._tidy_diag(c(HLR_chisq = hl$stat,
                                HLR_df = hl$df,
                                HLR_pvalue = hl$p.value)),
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
#'   \item{residuals}{Tibble with observed, fitted, Pearson residual.}
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
#' res = reg_poisson(y ~ x1, data = data.frame(y, x1))
#' res$coefficient
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
  obs_vec = data[[resp_col]]
  if (is.null(obs_vec)) obs_vec = as.numeric(stats::model.response(model$model))
  resid_df = tibble::tibble(
    observed         = obs_vec,
    fitted_val       = as.numeric(fitted(model)),
    pearson_residual = as.numeric(pearson_resid)
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
#'   \item{residuals}{Tibble of observed, fitted, Pearson residual.}
#' }
#'
#' @details
#' Unlike \code{reg_poisson()}, the negative binomial model handles
#' overdispersion through the estimated dispersion parameter \code{theta}
#' (fitted by \code{MASS::glm.nb()}): a small \code{theta} indicates strong
#' overdispersion, while \code{theta} growing large approaches the Poisson model.
#'
#' @examples
#' set.seed(42)
#' x1 = rnorm(100)
#' mu = exp(1 + 0.3*x1)
#' y = rnbinom(100, size = 1, mu = mu)
#' res = reg_negbin(y ~ x1, data = data.frame(y, x1))
#' res$coefficient
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
  obs_vec = data[[resp_col]]
  if (is.null(obs_vec)) obs_vec = as.numeric(stats::model.response(model$model))
  resid_df = tibble::tibble(
    observed         = obs_vec,
    fitted_val       = as.numeric(fitted(model)),
    pearson_residual = as.numeric(residuals(model, type = "pearson"))
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
      hl = tryCatch(
        .hosmer_lemeshow(model_result$residuals$observed,
                         model_result$residuals$fitted_val),
        error = function(e) list(stat = NA_real_, df = NA_integer_,
                                 p.value = NA_real_)
      )
      res$hosmer_lemeshow = tibble::tibble(
        chi_sq = round(hl$stat, 4), df = hl$df, p.value = round(hl$p.value, 6))
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
    sw = tryCatch(shapiro.test(resid_vec), error = function(e) NULL)
    res$shapiro_wilk = if (!is.null(sw)) {
      tibble::tibble(
        statistic = round(sw$statistic, 4),
        p.value   = round(sw$p.value, 6))
    } else {
      tibble::tibble(statistic = NA_real_, p.value = NA_real_)
    }

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
#' @param n_new Number of new observations to generate when \code{newdata}
#'   is not supplied. Predictor values are drawn from the empirical
#'   distributions of the input columns (normal approximation by mean/sd for
#'   numeric columns; levels sampled with observed proportions for
#'   factor columns).
#' @param level Confidence/interval level. Defaults to 0.95.
#' @param newdata Optional data frame of new predictor values. When supplied,
#'   predictions are made for exactly these rows and \code{n_new} is ignored.
#'
#' @return A named list with:
#' \describe{
#'   \item{predictions}{Tibble with predicted values.}
#'   \item{confidence_interval}{Tibble with lower and upper bounds of the confidence interval.}
#'   \item{prediction_interval}{Tibble with lower and upper bounds of the prediction interval (linear models only).}
#'   \item{newdata}{The data frame actually used for prediction.}
#' }
#' @export
reg_predict = function(model_result, n_new = 10, level = 0.95, newdata = NULL) {
  m = model_result$model

  if (inherits(m, "glm")) {
    fam = m$family
    if (fam$link == "logit") {
      return(.predict_logistic(m, model_result, n_new = n_new, level = level,
                               newdata = newdata))
    } else if (fam$link == "log") {
      return(.predict_count(m, model_result, n_new = n_new, level = level,
                            newdata = newdata))
    }
    stop("Unrecognized GLM family.", call. = FALSE)
  } else if (inherits(m, "lm")) {
    return(.predict_lm(m, model_result, n_new = n_new, level = level,
                       newdata = newdata))
  } else {
    stop("Unsupported model type.", call. = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Prediction helpers
# -----------------------------------------------------------------------------

.predict_lm = function(model, model_result, n_new, level, newdata = NULL) {
  newdf = .make_newdata(model, model_result, n_new, newdata = newdata)
  n_pred = nrow(newdf)
  pred_ci = predict(model, newdata = newdf, interval = "confidence", level = level)
  pred_pi = predict(model, newdata = newdf, interval = "prediction", level = level)

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_pred),
      .predict = as.numeric(pred_ci[, "fit"])),
    confidence_interval = tibble::tibble(
      lower = round(as.numeric(pred_ci[, "lwr"]), 4),
      upper = round(as.numeric(pred_ci[, "upr"]), 4)),
    prediction_interval = tibble::tibble(
      lower = round(as.numeric(pred_pi[, "lwr"]), 4),
      upper = round(as.numeric(pred_pi[, "upr"]), 4)),
    newdata = newdf
  )
}

.predict_logistic = function(model, model_result, n_new, level, newdata = NULL) {
  newdf = .make_newdata(model, model_result, n_new, newdata = newdata)
  n_pred = nrow(newdf)
  pred = predict(model, newdata = newdf, type = "link", se.fit = TRUE)
  z = qnorm(1 - (1 - level) / 2)
  fit_response = plogis(pred$fit)
  ci_lower = plogis(pred$fit - z * pred$se.fit)
  ci_upper = plogis(pred$fit + z * pred$se.fit)

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_pred),
      .predict = round(as.numeric(fit_response), 4)),
    confidence_interval = tibble::tibble(
      lower = round(ci_lower, 4),
      upper = round(ci_upper, 4)),
    newdata = newdf
  )
}

.predict_count = function(model, model_result, n_new, level, newdata = NULL) {
  newdf = .make_newdata(model, model_result, n_new, newdata = newdata)
  n_pred = nrow(newdf)
  pred = predict(model, newdata = newdf, type = "response", se.fit = TRUE)
  z = qnorm(1 - (1 - level) / 2)
  ci_lower = pmax(pred$fit - z * pred$se.fit, 0)
  ci_upper = pred$fit + z * pred$se.fit

  list(
    predictions = tibble::tibble(
      .new_id  = seq_len(n_pred),
      .predict = round(as.numeric(pred$fit), 4)),
    confidence_interval = tibble::tibble(
      lower = round(as.numeric(ci_lower), 4),
      upper = round(as.numeric(ci_upper), 4)),
    newdata = newdf
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
      title = paste(tryCatch(
        .reg_formula(model_result$formula)$response,
        error = function(e) "Model"), "Model Diagnostics"),
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
plot_reg_residuals = function(model_result) {
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
plot_reg_predict = function(model_result,
                            type = c("fitted", "predicted"),
                            show_ci = TRUE) {
  .prediction_plot_generic(model_result, type, show_ci = show_ci)
}
