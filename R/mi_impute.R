# ============================================================================
# Multiple imputation via a hand-written FCS (fully conditional
# specification) engine. Numeric columns use strict proper Bayesian
# linear-model draws (pmm / norm); binary and polytomous columns use
# parameter-perturbed logreg / polyreg (approximately proper, following
# the same spirit as mice).
# ============================================================================

# Number of donors used by the pmm matcher (mean of the 3 closest).
.FCS_PMM_DONORS = 3L

# Decide the default imputation method per column: numeric -> pmm;
# factor/character/logical with <= 2 levels -> logreg; > 2 -> polyreg.
.fcs_methods = function(data, col_nms) {
  vapply(col_nms, \(nm) {
    x = data[[nm]]
    if (is.numeric(x)) return("pmm")
    if (is.logical(x) || is.factor(x) || is.character(x)) {
      lv = if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
      return(if (length(lv) <= 2L) "logreg" else "polyreg")
    }
    stop("unsupported column type for imputation: ", nm, " (",
         class(x)[1], ").", call. = FALSE)
  }, character(1))
}

# Initial fill for the chain start: column mean (numeric) / mode
# (categorical). Character columns are coerced to factor in the working
# copy the chain iterates on.
.fcs_init = function(data, col_nms) {
  out = data
  for (nm in col_nms) {
    if (is.character(out[[nm]])) out[[nm]] = factor(out[[nm]])
  }
  for (nm in col_nms) {
    x = out[[nm]]
    ok = !is.na(x)
    if (all(ok)) next
    if (is.numeric(x)) {
      x[!ok] = mean(x[ok])
    } else {
      tab = table(x[ok])
      x[!ok] = names(tab)[which.max(tab)]
    }
    out[[nm]] = x
  }
  out
}

# Build design (features) and response for one target column. Features
# are the other columns that have no NA in df_cur (imputed columns are
# filled there; never-imputed NA columns are excluded so model.matrix
# cannot silently drop rows and misalign the indices). X carries an
# intercept as its first column.
.fcs_design = function(df_cur, data, nm) {
  na_idx = which(is.na(data[[nm]]))
  obs_idx = which(!is.na(data[[nm]]))
  others = setdiff(names(data), nm)
  others = others[!vapply(df_cur[others], anyNA, logical(1))]
  X = if (length(others) == 0L) {
    matrix(1, nrow(df_cur), 1)
  } else {
    stats::model.matrix(~ ., data = df_cur[, others, drop = FALSE])
  }
  list(X = X, y_obs = df_cur[[nm]][obs_idx], na_idx = na_idx,
       obs_idx = obs_idx)
}

# Bayesian linear-model draw with a non-informative prior: sigma2 ~
# scaled inverse chi-square on (n - p) df, then beta | sigma2 ~
# N(beta_hat, sigma2 (X'X)^-1). Ridge fallback for a singular
# cross-product.
.fcs_bayes_lm = function(X, y) {
  n = nrow(X)
  p = ncol(X)
  XTX_inv = tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XTX_inv)) {
    XTX_inv = solve(crossprod(X) + diag(1e-8, p))
  }
  beta_hat = as.numeric(XTX_inv %*% crossprod(X, y))
  resid = y - as.numeric(X %*% beta_hat)
  df_res = max(n - p, 1L)
  ss = sum(resid^2) / df_res
  sigma2 = ss * df_res / stats::rchisq(1, df_res)
  beta = beta_hat + as.numeric(
    t(chol(sigma2 * XTX_inv)) %*% matrix(stats::rnorm(p), p, 1))
  list(beta = beta, sigma2 = sigma2)
}

# Predicted-mean matching: fill each missing cell with the mean of the
# closest .FCS_PMM_DONORS observed values (deterministic donor variant).
.fcs_pmm = function(des) {
  pred_obs = as.numeric(des$X[des$obs_idx, , drop = FALSE] %*% des$beta)
  pred_na = as.numeric(des$X[des$na_idx, , drop = FALSE] %*% des$beta)
  y_obs = des$y_obs
  vapply(seq_along(pred_na), \(i) {
    d = abs(pred_obs - pred_na[i])
    nn = order(d)[seq_len(min(.FCS_PMM_DONORS, length(d)))]
    mean(y_obs[nn])
  }, numeric(1))
}

# One-time singularity warning via an environment (reference semantics):
# .warn_singular(rng_warn, "logreg", nm) warns only once per rng_warn.
.warn_singular = function(rng_warn, method, nm) {
  if (is.null(rng_warn$warned)) {
    warning("method '", method, "' for column '", nm,
            "': singular covariance; using unperturbed fit.",
            call. = FALSE)
    rng_warn$warned = TRUE
  }
}

# One FCS column update. Returns the refilled df_cur[[nm]].
.fcs_impute_col = function(df_cur, data, nm, method, rng_warn) {
  des = .fcs_design(df_cur, data, nm)
  if (length(des$na_idx) == 0L) return(df_cur[[nm]])

  if (method == "pmm" || method == "norm") {
    draw = .fcs_bayes_lm(des$X[des$obs_idx, , drop = FALSE], des$y_obs)
    if (method == "pmm") {
      des$beta = draw$beta
      fill = .fcs_pmm(des)
    } else {
      pred = as.numeric(des$X[des$na_idx, , drop = FALSE] %*% draw$beta)
      fill = pred + stats::rnorm(length(pred), sd = sqrt(draw$sigma2))
    }
    x = df_cur[[nm]]
    x[des$na_idx] = fill
    return(x)
  }

  if (method == "logreg") {
    y_fac = df_cur[[nm]]
    y_bin = unclass(y_fac) - 1L  # second level = 1
    Xo = des$X[des$obs_idx, , drop = FALSE]
    fit = tryCatch(
      suppressWarnings(
        stats::glm.fit(x = Xo, y = y_bin[des$obs_idx],
                       family = stats::binomial())
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || anyNA(fit$coefficients)) {
      return(df_cur[[nm]])  # keep last-iteration values
    }
    # information matrix from the final IRLS working weights
    V = tryCatch(
      solve(crossprod(Xo * sqrt(fit$weights))),
      error = function(e) NULL
    )
    beta = if (is.null(V)) {
      .warn_singular(rng_warn, "logreg", nm)
      fit$coefficients
    } else {
      as.numeric(fit$coefficients + t(chol(V)) %*%
        matrix(stats::rnorm(length(fit$coefficients)),
               length(fit$coefficients), 1))
    }
    lin = as.numeric(des$X[des$na_idx, , drop = FALSE] %*% beta)
    prob = stats::plogis(lin)
    fill = ifelse(stats::rbinom(length(lin), 1, prob) == 1L,
                  levels(y_fac)[2], levels(y_fac)[1])
    x = df_cur[[nm]]
    x[des$na_idx] = fill
    return(x)
  }

  # polyreg: multinomial logit + parameter perturbation
  y_fac = df_cur[[nm]]
  lvs = levels(y_fac)
  Xo = des$X[des$obs_idx, , drop = FALSE]
  dat_tr = data.frame(
    ..mi_y = y_fac[des$obs_idx],
    Xo[, -1, drop = FALSE],   # features without our intercept column
    check.names = FALSE
  )
  fit = tryCatch(
    nnet::multinom(..mi_y ~ ., data = dat_tr, trace = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || anyNA(as.numeric(stats::coef(fit)))) {
    return(df_cur[[nm]])
  }
  beta_hat = as.numeric(stats::coef(fit))  # (L-1) x p, column-major
  V = tryCatch(stats::vcov(fit), error = function(e) NULL)
  beta = if (is.null(V)) {
    .warn_singular(rng_warn, "polyreg", nm)
    beta_hat
  } else {
    tryCatch(
      beta_hat + as.numeric(t(chol(V + diag(1e-10, nrow(V)))) %*%
        matrix(stats::rnorm(nrow(V)), nrow(V), 1)),
      error = function(e) {
        .warn_singular(rng_warn, "polyreg", nm)
        beta_hat
      }
    )
  }
  Xna = des$X[des$na_idx, , drop = FALSE]
  B = matrix(beta, nrow = length(lvs) - 1L, ncol = ncol(Xna))
  eta = cbind(0, Xna %*% t(B))          # baseline level gets eta = 0
  mx = apply(eta, 1, max)               # numerically stable softmax
  prob = exp(eta - mx)
  prob = prob / rowSums(prob)
  idx = vapply(seq_len(nrow(prob)), \(i)
    sample.int(length(lvs), 1L, prob = prob[i, ]),
    integer(1))
  fill = lvs[idx]
  x = df_cur[[nm]]
  x[des$na_idx] = fill
  x
}
