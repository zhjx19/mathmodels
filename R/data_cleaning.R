# ============================================================================
# Data cleaning toolkit: missing-value diagnostics, imputation, and
# outlier handling. Tidy data frames in, tidy data frames out.
# ============================================================================

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

# Record changed cells (variable, row, old, new) as a tidy tibble.
# old/new are stored as characters so that mixed-type columns (numeric,
# character, factor) can be bound into one audit table.
.dc_changes = function(var, rows, old, new) {
  tibble::tibble(
    variable = var,
    row = as.integer(rows),
    old = as.character(old),
    new = as.character(new)
  )
}

# Attach changes (and optional details) attributes to a returned data frame.
.dc_with_changes = function(df, changes, details = NULL) {
  df = structure(df, changes = changes)
  if (!is.null(details)) df = structure(df, details = details)
  df
}

# ----------------------------------------------------------------------------
# Missing-value diagnostics
# ----------------------------------------------------------------------------

#' Missing Value Summary
#'
#' Summarizes missing values per column of a tidy data frame: the number and
#' percentage of missing entries. The row dimension is covered by
#' [plot_na_heatmap()].
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Columns to summarize;
#'   defaults to all columns.
#'
#' @return A tibble: `variable`, `n_missing`, `pct_missing` (original column
#'   order).
#'
#' @examples
#' d = data.frame(x = c(1, NA, 3), y = c(NA, NA, 2))
#' na_summary(d)
#' @export
na_summary = function(data, .cols = dplyr::everything()) {
  .mv_frame(data)
  sel = dplyr::select(data, {{ .cols }})
  out = tibble::tibble(
    variable = names(sel),
    n_missing = unname(vapply(sel, \(x) sum(is.na(x)), integer(1))),
    pct_missing = unname(vapply(sel, \(x) mean(is.na(x)), numeric(1)))
  )
  out
}

# EM estimate of the mean vector and covariance matrix under the MCAR
# assumption (multivariate normal; Little & Rubin). Each missingness
# pattern fills its missing cells with the conditional expectations given
# the observed cells and accumulates the conditional covariance terms.
.em_norm_mle = function(X, maxit = 200, tol = 1e-8) {
  n = nrow(X)
  p = ncol(X)
  mu = colMeans(X, na.rm = TRUE)
  if (anyNA(mu)) stop("a column is fully missing.", call. = FALSE)
  sigma = stats::cov(X, use = "pairwise.complete.obs")
  if (anyNA(sigma)) stop("covariance cannot be estimated.", call. = FALSE)
  miss = is.na(X)
  pat_keys = apply(miss, 1, paste, collapse = "")
  groups = split(seq_len(n), pat_keys)
  patterns = lapply(groups, function(gi) miss[gi[1], ])

  for (it in seq_len(maxit)) {
    mu_old = mu
    sigma_old = sigma
    sum_x = matrix(0, p, p)
    sum_vec = numeric(p)
    for (ks in seq_along(groups)) {
      gi = groups[[ks]]
      mp = patterns[[ks]]
      m_cols = which(mp)
      o_cols = which(!mp)
      if (length(m_cols) == 0) {
        Xg = X[gi, , drop = FALSE]
        sum_vec = sum_vec + colSums(Xg)
        sum_x = sum_x + crossprod(Xg)
        next
      }
      S_oo = sigma[o_cols, o_cols, drop = FALSE]
      S_mo = sigma[m_cols, o_cols, drop = FALSE]
      # E step: conditional means and covariance of the missing block
      b = S_mo %*% solve(S_oo)
      Xg = X[gi, , drop = FALSE]
      Xf = Xg
      Xf[, m_cols] = sweep(Xg[, o_cols, drop = FALSE], 2, mu[o_cols]) %*%
        t(b) + matrix(mu[m_cols], nrow = length(gi), ncol = length(m_cols),
          byrow = TRUE)
      var_m = sigma[m_cols, m_cols, drop = FALSE] - b %*% t(S_mo)
      sum_vec = sum_vec + colSums(Xf)
      # accumulate E[x x']: filled block plus the conditional variance
      # correction for the missing block
      sum_x = sum_x + crossprod(Xf) + length(gi) *
        .em_block_var(var_m, m_cols, p)
    }
    mu = sum_vec / n
    sigma = sum_x / n - outer(mu, mu)
    # ensure symmetry after numerical noise
    sigma = (sigma + t(sigma)) / 2
    if (max(abs(mu - mu_old)) < tol &&
        max(abs(sigma - sigma_old)) < tol * max(1, max(abs(sigma)))) {
      break
    }
  }
  list(mu = mu, sigma = sigma, iterations = it)
}

# Expand a conditional-covariance block into a full p x p matrix (zeros
# outside the missing columns).
.em_block_var = function(var_m, m_cols, p) {
  out = matrix(0, p, p)
  out[m_cols, m_cols] = var_m
  out
}

#' Little's MCAR Test
#'
#' Performs Little's (1988) test of Missing Completely At Random (the
#' implementation follows `naniar::mcar_test()` / Hawkins' refinement):
#' the mean vector and covariance matrix are first estimated under the
#' MCAR assumption by maximum likelihood (EM algorithm), then each
#' missingness pattern's mean (on its observed columns) is compared with
#' the overall mean via a Mahalanobis distance weighted by the pattern
#' size. The statistic is referred to a chi-square distribution with
#' `df = sum(observed columns per pattern) - p`.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric columns to
#'   test; defaults to all columns.
#' @param maxit Maximum EM iterations.
#' @param tol EM convergence tolerance.
#'
#' @return A one-row tibble: `chi2`, `df`, `p_value`, `n_patterns`, `n_used`,
#'   `p_missing_cells`. Attributes: `patterns` — a tibble with one row per
#'   missingness pattern (missing columns, number of rows).
#'
#' @details
#' The EM estimate requires an invertible covariance matrix; fully missing
#' columns are dropped with a warning before fitting. With no missing
#' values (or a degenerate pattern structure) the test is not applicable
#' and an error is raised.
#'
#' @examples
#' set.seed(1)
#' d = data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
#' d$y[1:10] = NA
#' na_mcar_test(d)
#' @export
na_mcar_test = function(data, .cols = dplyr::everything(), maxit = 200,
                        tol = 1e-8) {
  .mv_frame(data)
  sel = dplyr::select(data, {{ .cols }})
  bad = names(sel)[!vapply(sel, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(".cols must select numeric column(s). Non-numeric: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  X = as.matrix(sel)
  n_na_col = colSums(is.na(X))
  if (all(n_na_col == 0)) {
    stop("data has no missing values; the MCAR test is not applicable.",
         call. = FALSE)
  }
  if (any(n_na_col == nrow(X))) {
    fully = names(n_na_col)[n_na_col == nrow(X)]
    warning("fully missing column(s) dropped: ",
            paste(fully, collapse = ", "), ".", call. = FALSE)
    X = X[, n_na_col < nrow(X), drop = FALSE]
  }
  p = ncol(X)
  em = .em_norm_mle(X, maxit = maxit, tol = tol)
  mu = em$mu
  sigma_inv = tryCatch(solve(em$sigma), error = function(e) NULL)
  if (is.null(sigma_inv)) {
    stop("the EM-estimated covariance matrix is singular; drop redundant ",
         "or highly-missing variables.", call. = FALSE)
  }

  miss = is.na(X)
  pat_keys = apply(miss, 1, paste, collapse = "")
  groups = split(seq_len(nrow(X)), pat_keys)
  chi2 = 0
  df = 0L
  pat_tbl = tibble::tibble(
    pattern = character(), n_rows = integer(), n_missing_cols = integer()
  )
  for (g in groups) {
    mp = miss[g[1], ]
    obs_cols = which(!mp)
    # Little (1988): compare each pattern's mean (on its observed columns)
    # with the overall (EM) mean vector, weighted by the pattern size
    mu_g = colMeans(X[g, obs_cols, drop = FALSE])
    dif = mu_g - mu[obs_cols]
    d2_g = as.numeric(t(dif) %*% sigma_inv[obs_cols, obs_cols,
      drop = FALSE] %*% dif)
    chi2 = chi2 + length(g) * d2_g
    df = df + length(obs_cols)
    pat_tbl = tibble::add_row(
      pat_tbl,
      pattern = paste(names(mp)[mp], collapse = ", "),
      n_rows = length(g),
      n_missing_cols = sum(mp)
    )
  }
  df = df - p
  p_value = if (df > 0) stats::pchisq(chi2, df = df, lower.tail = FALSE)
    else NA_real_
  out = tibble::tibble(
    chi2 = chi2,
    df = df,
    p_value = p_value,
    n_patterns = length(groups),
    n_used = nrow(X),
    p_missing_cells = mean(is.na(X))
  )
  structure(out, patterns = pat_tbl)
}

# ----------------------------------------------------------------------------
# Missing-value visualization
# ----------------------------------------------------------------------------

#' Missingness Heatmap
#'
#' Draws a row-by-column missingness heatmap: each cell shows whether the
#' entry is observed (grey) or missing (red). Rows keep their original order
#' so the plot stays aligned with the data.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Columns to display.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d = data.frame(x = c(1, NA, 3, 4), y = c(NA, 2, NA, 5))
#' plot_na_heatmap(d)
#' @export
plot_na_heatmap = function(data, .cols) {
  .mv_frame(data)
  sel = dplyr::select(data, {{ .cols }})
  m = is.na(as.matrix(sel))
  long = tibble::tibble(
    row = rep(seq_len(nrow(sel)), times = ncol(sel)),
    variable = factor(rep(names(sel), each = nrow(sel)),
      levels = rev(names(sel))),
    missing = as.vector(m)
  )
  ggplot2::ggplot(long,
      ggplot2::aes(x = .data$row, y = .data$variable, fill = .data$missing)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "grey90", `TRUE` = "#d62728"),
      labels = c(Observed = "Observed", Missing = "Missing")) +
    ggplot2::labs(x = "Row", y = NULL, fill = NULL) +
    ggplot2::theme_minimal()
}

#' Missingness Bar Plot
#'
#' Draws a bar plot of the missing-value percentage per column.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Columns to display.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d = data.frame(x = c(1, NA, 3, 4), y = c(NA, 2, NA, 5))
#' plot_na_bar(d)
#' @export
plot_na_bar = function(data, .cols) {
  summ = na_summary(data, {{ .cols }})
  summ$label = sprintf("%.1f%%", 100 * summ$pct_missing)
  ggplot2::ggplot(summ,
      ggplot2::aes(x = .data$variable, y = .data$pct_missing)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), vjust = -0.3) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
      labels = \(p) sprintf("%d%%", round(100 * p))) +
    ggplot2::labs(x = NULL, y = "Missing percentage") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal()
}

# ----------------------------------------------------------------------------
# Simple imputation
# ----------------------------------------------------------------------------

# Fill one column (vector xi, with NAs) using a simple method; returns the
# filled vector or NULL when nothing to fill.
.fill_simple = function(xi, method, value) {
  ok = !is.na(xi)
  n_ok = sum(ok)
  fill = switch(method,
    mean = if (n_ok) mean(xi[ok]) else NA,
    median = if (n_ok) stats::median(xi[ok]) else NA,
    mode = {
      tab = table(xi[ok])
      if (length(tab) == 0) NA else {
        m = names(tab)[which.max(tab)]
        if (is.numeric(xi)) as.numeric(m) else m
      }
    },
    constant = value
  )
  if (length(fill) != 1 || is.na(fill)) {
    warning("nothing to impute from for method = '", method,
      "'; column left as is.", call. = FALSE)
    return(NULL)
  }
  xi[is.na(xi)] = fill
  xi
}

#' Simple Missing Value Imputation
#'
#' Fills missing values in the selected columns with a simple statistic:
#' group mean, median, or mode (with `.by` slice grouping), row-order linear
#' or cubic-spline interpolation (reusing `interp_linear()` /
#' `interp_spline()`), or a constant. Returns a data frame of the same shape;
#' the changed cells are recorded in `attr(result, "changes")`.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Columns to impute.
#' @param method `"mean"`, `"median"`, `"mode"`, `"linear"`, `"spline"`, or
#'   `"constant"`.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns;
#'   the statistic is computed within each slice. Supported for `"mean"`,
#'   `"median"`, and `"mode"`.
#' @param value Constant fill value (required for `method = "constant"`).
#'
#' @return The data frame with missing values filled, with attribute
#'   `changes` (tibble: `variable`, `row`, `old`, `new`).
#'
#' @details
#' `"linear"` / `"spline"` interpolate over row order, i.e. they are meant
#' for ordered (time-series-like) data. Missing entries at the two ends
#' cannot be interpolated and remain `NA`. Fully missing columns stay `NA`
#' with a warning. `"mode"` works for any column type; `"mean"` and
#' `"median"` require numeric columns.
#'
#' @examples
#' d = data.frame(
#'   g = rep(c("a", "b"), each = 3),
#'   x = c(1, NA, 3, 10, NA, 12)
#' )
#' impute(d, .cols = x, method = "mean", .by = g)
#' impute(d, .cols = x, method = "constant", value = 0)
#' @export
impute = function(data, .cols,
                  method = c("mean", "median", "mode", "linear", "spline",
                             "constant"),
                  .by = NULL, value = NULL) {
  .mv_frame(data)
  method = match.arg(method)
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  col_nms = names(cols_df)
  if (method == "constant" && is.null(value)) {
    stop("method = 'constant' requires `value`.", call. = FALSE)
  }
  if (method %in% c("mean", "median", "linear", "spline")) {
    bad = col_nms[!vapply(cols_df, is.numeric, logical(1))]
    if (length(bad) > 0) {
      stop("method = '", method, "' requires numeric column(s). ",
           "Non-numeric: ", paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (method %in% c("linear", "spline") && ncol(slice_df) > 0) {
    stop("method = '", method, "' does not support .by grouping.",
         call. = FALSE)
  }

  sl = .slice_index(slice_df)
  changes_list = vector("list", length(sl$idx) * length(col_nms))
  filled_cols = vector("list", length(col_nms))
  pos = 1L
  for (nm in col_nms) {
    x = cols_df[[nm]]
    for (k in seq_along(sl$idx)) {
      i = sl$idx[[k]]
      xi = x[i]
      na_idx = which(is.na(xi))
      if (length(na_idx) == 0) next
      if (method %in% c("linear", "spline")) {
        ok = which(!is.na(xi))
        min_ok = if (method == "linear") 2L else 4L
        if (length(ok) < min_ok) {
          warning("column '", nm, "': too few observed values for ",
                  method, " interpolation; column left as is.",
                  call. = FALSE)
          next
        }
        res = if (method == "linear") {
          interp_linear(ok, xi[ok], seq_along(xi))
        } else {
          interp_spline(ok, xi[ok], seq_along(xi))
        }
        fitted = as.numeric(res$y)
        new_vals = fitted[na_idx]
        xi[na_idx] = new_vals
      } else {
        xi = .fill_simple(xi, method, value)
        if (is.null(xi)) next
        new_vals = xi[na_idx]
      }
      old_vals = x[i][na_idx]
      changes_list[[pos]] = .dc_changes(
        var = nm, rows = i[na_idx], old = old_vals, new = new_vals)
      pos = pos + 1L
      x[i] = xi
    }
    filled_cols[[nm]] = x
  }
  changes = dplyr::bind_rows(changes_list)
  out = data
  for (nm in col_nms) out[[nm]] = filled_cols[[nm]]
  .dc_with_changes(out, changes)
}

# ----------------------------------------------------------------------------
# Model-based imputation (unified framework: knn / lm / tree)
# ----------------------------------------------------------------------------

#' Model-Based Missing Value Imputation
#'
#' Unified imputation framework on top of three models: k-nearest-neighbors
#' (`knn`, hand-written), linear regression (`lm`, `stats::lm`), and
#' regression trees (`tree`, `rpart::rpart`). For every numeric column with
#' missing values, a model is trained on the complete rows (target column and
#' all other numeric columns as features) and the missing cells are replaced
#' by the predictions. Returns a data frame of the same shape with the
#' changed cells recorded in `attr(result, "changes")` and per-column details
#' in `attr(result, "details")`.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric target
#'   columns to impute. Features are all other numeric columns of `data`.
#' @param model `"knn"`, `"lm"`, or `"tree"`.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns;
#'   models are trained separately within each slice.
#' @param k Number of neighbors for `model = "knn"`.
#' @param ... Additional arguments passed to `rpart::rpart()` for
#'   `model = "tree"` (e.g. `control = rpart::rpart.control(cp = 0.01)`).
#'
#' @return The data frame with missing values filled, with attributes
#'   `changes` (tibble: `variable`, `row`, `old`, `new`) and `details`
#'   (tibble: `variable`, `model`, `n_train`, `n_imputed`).
#'
#' @details
#' Rows whose features still contain missing values cannot be scored by
#' `"lm"` / `"tree"` predictions and keep their `NA` (a joint iterative
#' scheme is planned for a later version). `"knn"` handles rows with
#' partially missing features by computing distances on the features
#' observed for the target row; rows with all features missing keep `NA`.
#'
#' @examples
#' set.seed(1)
#' d = data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
#' d$x[c(3, 10)] = NA
#' impute_model(d, .cols = x, model = "knn", k = 5)
#' impute_model(d, .cols = x, model = "lm")
#' @export
impute_model = function(data, .cols,
                        model = c("knn", "lm", "tree"),
                        .by = NULL, k = 5, ...) {
  .mv_frame(data)
  model = match.arg(model)
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  col_nms = names(cols_df)
  bad = col_nms[!vapply(cols_df, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(".cols must select numeric column(s). Non-numeric: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  if (model == "knn") {
    if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
      stop("k must be a positive integer.", call. = FALSE)
    }
  }
  # features: all other numeric columns of data
  feat_nms = setdiff(
    names(data)[vapply(data, is.numeric, logical(1))],
    col_nms
  )
  if (length(feat_nms) == 0) {
    stop("no numeric feature columns available outside .cols.",
         call. = FALSE)
  }

  sl = .slice_index(slice_df)
  changes_list = vector("list", length(sl$idx) * length(col_nms))
  details_list = vector("list", length(sl$idx) * length(col_nms))
  pos = 1L
  out = data
  for (k_slice in seq_along(sl$idx)) {
    i_slice = sl$idx[[k_slice]]
    for (nm in col_nms) {
      i = i_slice
      y_full = data[[nm]]
      train_idx = i[stats::complete.cases(data[i, c(nm, feat_nms),
        drop = FALSE])]
      if (length(train_idx) < 2) {
        warning("column '", nm, "': fewer than 2 complete rows; skipped.",
                call. = FALSE)
        next
      }
      y = y_full[train_idx]
      X_train = as.matrix(data[train_idx, feat_nms, drop = FALSE])
      na_idx = i[which(is.na(y_full[i]))]
      if (length(na_idx) == 0) next
      pred = rep(NA_real_, length(na_idx))
      if (model == "knn") {
        X_na = as.matrix(data[na_idx, feat_nms, drop = FALSE])
        kk = min(k, length(train_idx))
        for (r in seq_along(na_idx)) {
          feats = which(!is.na(X_na[r, ]))
          if (length(feats) == 0) next  # keep NA
          dist = sqrt(rowSums((X_train[, feats, drop = FALSE] -
            matrix(X_na[r, feats], nrow(X_train), length(feats),
                   byrow = TRUE))^2))
          nn = order(dist)[seq_len(kk)]
          pred[r] = mean(y[nn])
        }
      } else {
        train_df = data.frame(y = y, X_train)
        names(train_df) = c(".__y", feat_nms)
        f = stats::as.formula(paste0(".__y ~ ",
          paste(paste0("`", feat_nms, "`"), collapse = " + ")))
        fit = if (model == "lm") {
          stats::lm(f, data = train_df)
        } else {
          rpart::rpart(f, data = train_df, ...)
        }
        na_df = data[na_idx, feat_nms, drop = FALSE]
        names(na_df) = feat_nms
        pred = tryCatch(
          as.numeric(stats::predict(fit, newdata = na_df)),
          error = function(e) rep(NA_real_, nrow(na_df)))
      }
      old_vals = y_full[na_idx]
      ok_pred = !is.na(pred)
      if (any(ok_pred)) {
        out[na_idx[ok_pred], nm] = pred[ok_pred]
        changes_list[[pos]] = .dc_changes(
          var = nm, rows = na_idx[ok_pred],
          old = old_vals[ok_pred], new = pred[ok_pred])
      }
      details_list[[pos]] = tibble::tibble(
        variable = nm,
        model = model,
        n_train = length(train_idx),
        n_imputed = sum(ok_pred)
      )
      pos = pos + 1L
    }
  }
  changes = dplyr::bind_rows(changes_list)
  details = dplyr::bind_rows(details_list)
  .dc_with_changes(out, changes, details)
}

# ----------------------------------------------------------------------------
# Outliers: detect / replace with NA / cap
# ----------------------------------------------------------------------------

# Compute two-sided bounds for one column under a detection method.
.outlier_bounds = function(x, method, k) {
  switch(method,
    iqr = {
      q = stats::quantile(x, c(0.25, 0.75), names = FALSE, na.rm = TRUE)
      iqr = q[2] - q[1]
      if (iqr == 0) return(NULL)
      c(q[1] - k * iqr, q[2] + k * iqr)
    },
    zscore = {
      m = mean(x); s = stats::sd(x)
      if (s == 0 || is.na(s)) return(NULL)
      c(m - k * s, m + k * s)
    },
    mad = {
      med = stats::median(x)
      s = 1.4826 * stats::mad(x)
      if (s == 0 || is.na(s)) return(NULL)
      c(med - k * s, med + k * s)
    }
  )
}

#' Outlier Detection and Handling
#'
#' Flags outliers in the selected numeric columns with a two-sided bound rule
#' and, optionally, acts on them. Detection methods: Tukey's IQR fence
#' (`iqr`), z-score (`zscore`, `k` standard deviations around the mean), and
#' the Hampel filter (`mad`, `k` scaled-MAD distances around the median).
#' By default the data is left untouched (`to = "detect"`): the detection
#' detail table is returned. With `to = "na"` detected values become `NA`
#' (ready for `impute()`); with `to = "cap"` they are capped at the fence
#' values.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric columns.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns;
#'   bounds are computed within each slice.
#' @param method `"iqr"`, `"zscore"`, or `"mad"`.
#' @param k Fence multiplier. `NULL` (default) picks the conventional value
#'   per method: 1.5 for `iqr`, 3 for `zscore`, 3 for `mad`.
#' @param to `"detect"` (default; data untouched, detection table returned),
#'   `"na"` (detected values set to `NA`), or `"cap"` (capped at the bounds).
#'
#' @return
#' `to = "detect"`: a tibble with one row per detected outlier (`variable`,
#' `.by` identifiers, `row`, `value`, `bound_low`, `bound_high`), with
#' attribute `summary` (per-column counts). `to = "na"` / `"cap"`: the data
#' frame with attribute `outliers` (the same detail table) and `summary`.
#'
#' @details
#' Constant columns (zero IQR / sd / MAD) cannot be scored and are skipped
#' with a warning.
#'
#' @examples
#' d = data.frame(g = rep(c("a", "b"), each = 20),
#'                x = c(rnorm(20), rnorm(19, 10), 100))
#' outlier(d, .cols = x, .by = g)
#' outlier(d, .cols = x, .by = g, to = "na")
#' outlier(d, .cols = x, .by = g, to = "cap")
#' @export
outlier = function(data, .cols, .by = NULL,
                   method = c("iqr", "zscore", "mad"),
                   k = NULL, to = c("detect", "na", "cap")) {
  .mv_frame(data)
  method = match.arg(method)
  to = match.arg(to)
  if (is.null(k)) k = switch(method, iqr = 1.5, zscore = 3, mad = 3)
  if (!is.numeric(k) || length(k) != 1 || k <= 0) {
    stop("k must be a positive number.", call. = FALSE)
  }
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  col_nms = names(cols_df)
  bad = col_nms[!vapply(cols_df, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(".cols must select numeric column(s). Non-numeric: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }

  sl = .slice_index(slice_df)
  detail_list = vector("list", length(sl$idx) * length(col_nms))
  pos = 1L
  out = data
  skipped = character(0)
  for (ks in seq_along(sl$idx)) {
    i = sl$idx[[ks]]
    for (nm in col_nms) {
      x = data[[nm]][i]
      ok = !is.na(x)
      bounds = .outlier_bounds(x[ok], method, k)
      if (is.null(bounds)) {
        skipped = c(skipped, nm)
        next
      }
      flag = ok & (x < bounds[1] | x > bounds[2])
      idx_hit = i[flag]
      if (length(idx_hit) > 0) {
        detail_list[[pos]] = tibble::tibble(
          variable = nm,
          row = idx_hit,
          value = data[[nm]][idx_hit],
          bound_low = bounds[1],
          bound_high = bounds[2]
        )
        pos = pos + 1L
        if (to == "na") {
          out[idx_hit, nm] = NA
        } else if (to == "cap") {
          out[idx_hit, nm] = pmin(pmax(data[[nm]][idx_hit], bounds[1]),
            bounds[2])
        }
      }
    }
  }
  if (length(skipped) > 0) {
    warning("constant column(s) skipped (zero spread): ",
            paste(unique(skipped), collapse = ", "), ".", call. = FALSE)
  }
  detail = dplyr::bind_rows(detail_list)
  # slice keys: map each detected row to its slice key
  if (ncol(slice_df) > 0 && nrow(detail) > 0) {
    slice_of_row = integer(nrow(data))
    for (ks in seq_along(sl$idx)) slice_of_row[sl$idx[[ks]]] = ks
    keys = sl$keys[slice_of_row[detail$row], , drop = FALSE]
    detail = dplyr::bind_cols(keys, detail)
    detail = detail[, c(names(slice_df),
      setdiff(names(detail), names(slice_df))), drop = FALSE]
  }
  summary = if (nrow(detail) > 0) {
    detail |>
      dplyr::summarise(
        n_outliers = dplyr::n(),
        .by = c(dplyr::any_of(names(slice_df)), dplyr::all_of("variable")))
  } else {
    tibble::tibble(variable = character(), n_outliers = integer())
  }
  summary$method = method
  if (to == "detect") {
    structure(detail, summary = summary)
  } else {
    structure(out, outliers = detail, summary = summary)
  }
}

# ----------------------------------------------------------------------------
# Pivot table
# ----------------------------------------------------------------------------

#' Pivot Table
#'
#' Builds a pivot table: aggregates `values` by the cross of `index` rows and
#' `columns`, then spreads the aggregates into a wide table (one column per
#' column level x value).
#'
#' @param data A data frame.
#' @param index <[`tidy-select`][tidyselect::language]> Row identifier
#'   column(s).
#' @param columns <[`tidy-select`][tidyselect::language]> Column factor(s)
#'   to spread over.
#' @param values <[`tidy-select`][tidyselect::language]> Numeric value
#'   column(s) to aggregate.
#' @param aggfunc Aggregation function (default `mean`).
#' @param fill_value Value used for empty cells (default `NA`).
#' @param names_sep Separator between column names and value names when
#'   several `values` are spread.
#'
#' @return A tibble in wide form.
#'
#' @examples
#' pivot_table(ToothGrowth, index = supp, columns = dose, values = len)
#' pivot_table(ToothGrowth, index = supp, columns = dose, values = len,
#'   aggfunc = length)
#' @export
pivot_table = function(data, index, columns, values, aggfunc = mean,
                       fill_value = NA, names_sep = "_") {
  .mv_frame(data)
  nm_cols = names(dplyr::select(data, {{ columns }}))
  names_glue = paste0("{.value}", names_sep,
    "{", paste(nm_cols, collapse = names_sep), "}")
  data |>
    dplyr::summarise(
      dplyr::across({{ values }}, aggfunc),
      .by = c({{ index }}, {{ columns }})) |>
    tidyr::pivot_wider(
      id_cols = {{ index }},
      names_from = {{ columns }},
      values_from = {{ values }},
      names_glue = names_glue,
      values_fill = fill_value
    )
}
