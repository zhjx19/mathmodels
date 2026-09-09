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
# (categorical). Character and logical columns are coerced to factor
# in the working copy the chain iterates on (logicals must become
# factors before the logreg path's unclass() indexing).
.fcs_init = function(data, col_nms) {
  out = data
  for (nm in col_nms) {
    if (is.character(out[[nm]]) || is.logical(out[[nm]])) {
      out[[nm]] = factor(out[[nm]])
    }
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
  if (!(sigma2 > 0)) {
    # degenerate fit (n <= p or zero residuals): sigma2 -> 0 would make
    # chol() see a zero matrix; fall back to the unperturbed OLS draw
    return(list(beta = beta_hat, sigma2 = 0))
  }
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
  # coef(fit) is (L-1) x p with rows = levels; unwinding it directly
  # (as.numeric(coef)) is PREDICTOR-major, but stats::vcov(fit) orders
  # coefficients LEVEL-major. Stack beta_hat level-major (via t()) so
  # each coefficient's perturbation draws from its own variance slot.
  beta_hat = as.numeric(t(stats::coef(fit)))  # level-major, matches vcov
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
  # beta is level-major: matrix(beta, nrow = p) fills predictor-major
  # (p x (L-1)) from it, so t() gives (L-1) x p with rows = levels in
  # vcov order; cbind(0, .) puts the baseline level's eta at 0.
  B = t(matrix(beta, nrow = ncol(Xna)))
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

# ----------------------------------------------------------------------------
# Exported interface
# ----------------------------------------------------------------------------

#' Multiple Imputation via FCS
#'
#' Generates `m` completed data sets by fully conditional specification
#' (chained equations). Numeric columns use proper Bayesian linear-model
#' draws with predicted-mean matching (`"pmm"`, default) or normal noise
#' (`"norm"`); binary and polytomous columns use parameter-perturbed
#' logistic (`"logreg"`) / multinomial (`"polyreg"`) models — approximately
#' proper, following the same spirit as `mice`.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Columns to impute;
#'   defaults to all columns. Columns without missing values are left as is.
#' @param m Number of imputations (at least 2).
#' @param maxit Number of chain iterations per imputation.
#' @param seed Random seed for reproducibility.
#' @param methods Optional named vector/list overriding the automatic
#'   method choice per column, e.g. `methods = c(x = "norm")`.
#'
#' @return An object of class `mi_impute`: a list of `m` completed data
#'   frames. Attributes: `methods`, `iterations`, `trace` (tibble with
#'   imputation / iteration / variable / mean / sd of the imputed cells;
#'   mean/sd are NA for categorical columns), `seed`, `n_imputed`.
#'
#' @details
#' Logical and character columns are treated as categorical and come back
#' as factors. Convergence can be inspected via the `trace` attribute
#' (plotting is planned for a later version).
#'
#' @examples
#' set.seed(1)
#' d = data.frame(x = rnorm(30), z = rnorm(30),
#'                g = factor(sample(c("a", "b"), 30, replace = TRUE)))
#' d$x[c(2, 8)] = NA
#' d$g[5] = NA
#' imps = impute_multiple(d, m = 3, maxit = 3, seed = 11)
#' imps
#' @export
impute_multiple = function(data, .cols = dplyr::everything(),
                           m = 5, maxit = 5, seed = NULL, methods = NULL) {
  .mv_frame(data)
  if (!is.numeric(m) || length(m) != 1L || m < 2 || m != floor(m)) {
    stop("m must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1 ||
      maxit != floor(maxit)) {
    stop("maxit must be an integer >= 1.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) {
      stop("seed must be a single number.", call. = FALSE)
    }
    set.seed(seed)
  }

  cols_df = tryCatch(
    dplyr::select(data, {{ .cols }}),
    error = function(e) NULL
  )
  if (is.null(cols_df)) {
    stop("`.cols` selected no existing column(s).", call. = FALSE)
  }
  col_nms = names(cols_df)
  # only columns that actually need imputation
  n_na_col = vapply(data[col_nms], \(x) sum(is.na(x)), integer(1))
  col_nms = col_nms[n_na_col > 0]
  # fully missing columns have no observed values to model from: leave
  # them NA with a warning (same convention as impute())
  full_miss = col_nms[n_na_col[col_nms] == nrow(data)]
  if (length(full_miss) > 0) {
    warning("fully missing column(s) left as NA: ",
            paste(full_miss, collapse = ", "), ".", call. = FALSE)
    col_nms = setdiff(col_nms, full_miss)
  }
  if (length(col_nms) == 0L) {
    res = unclass(rep(list(data), m))
    class(res) = c("mi_impute", "list")
    attr(res, "methods") = setNames(character(), character())
    attr(res, "iterations") = 0L
    attr(res, "trace") = tibble::tibble(
      imputation = integer(), iteration = integer(),
      variable = character(), mean = numeric(), sd = numeric()
    )
    attr(res, "seed") = seed
    attr(res, "n_imputed") = setNames(integer(), character())
    return(res)
  }

  # resolve methods: automatic per type, then user overrides
  auto = .fcs_methods(data, col_nms)
  if (!is.null(methods)) {
    if (!is.list(methods)) methods = as.list(methods)
    unknown = setdiff(unlist(methods, use.names = FALSE),
                      c("pmm", "norm", "logreg", "polyreg"))
    if (length(unknown) > 0) {
      stop("unknown imputation method(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    for (nm in names(methods)) {
      if (!(nm %in% col_nms)) {
        stop("methods names a column that is not imputed: ", nm,
             ".", call. = FALSE)
      }
      auto[[nm]] = methods[[nm]]
    }
  }

  out_sets = vector("list", m)
  trace_rows = vector("list", m * maxit * length(col_nms))
  pos = 1L
  for (mi in seq_len(m)) {
    rng_warn = new.env(parent = emptyenv())
    df_cur = .fcs_init(data, col_nms)
    for (it in seq_len(maxit)) {
      for (nm in col_nms) {
        df_cur[[nm]] = .fcs_impute_col(df_cur, data, nm,
                                       method = auto[[nm]], rng_warn)
      }
      for (nm in col_nms) {
        vals = df_cur[[nm]][is.na(data[[nm]])]
        trace_rows[[pos]] = tibble::tibble(
          imputation = mi, iteration = it, variable = nm,
          mean = if (is.numeric(vals)) mean(vals) else NA_real_,
          sd = if (is.numeric(vals)) stats::sd(vals) else NA_real_
        )
        pos = pos + 1L
      }
    }
    out_sets[[mi]] = df_cur
  }
  trace = dplyr::bind_rows(trace_rows)

  res = out_sets
  class(res) = c("mi_impute", "list")
  attr(res, "methods") = auto
  attr(res, "iterations") = maxit
  attr(res, "trace") = trace
  attr(res, "seed") = seed
  attr(res, "n_imputed") = n_na_col[col_nms]
  res
}

#' @export
print.mi_impute = function(x, ...) {
  m = length(x)
  mth = attr(x, "methods")
  cat("Multiple imputation (mi_impute):", m, "imputation(s) of",
      length(mth), "column(s)\n")
  if (length(mth) > 0) {
    cat("  methods:",
        paste(sprintf("%s=%s", names(mth), mth), collapse = ", "), "\n")
  }
  ni = attr(x, "n_imputed")
  if (length(ni) > 0) {
    cat("  cells imputed:",
        paste(sprintf("%s=%d", names(ni), ni), collapse = ", "), "\n")
  }
  cat("  iterations:", attr(x, "iterations"), "\n")
  invisible(x)
}
