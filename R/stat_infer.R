# ============================================================================
# Statistical inference toolkit: tidy-in, tidy-out hypothesis tests,
# bootstrap and permutation frameworks (pure tidy data frames in and out)
# ============================================================================

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

# Validate that the input is a data frame (same contract as .mv_frame).
.stat_frame = function(data) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  invisible(TRUE)
}

# Significance stars from p-values.
.sig_stars = function(p) {
  out = rep("", length(p))
  out[!is.na(p) & p < 0.1] = "."
  out[!is.na(p) & p < 0.05] = "*"
  out[!is.na(p) & p < 0.01] = "**"
  out[!is.na(p) & p < 0.001] = "***"
  out
}

# Validate alternative / conf_level / p_adjust arguments shared across tests.
.check_alt = function(alternative) {
  match.arg(alternative, c("two.sided", "greater", "less"))
}

.check_conf = function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single value in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

.check_padjust = function(p_adjust) {
  if (length(p_adjust) != 1 || !p_adjust %in% stats::p.adjust.methods) {
    stop("p_adjust must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

# Split row indices by the slice columns (.by). Returns list(idx, keys):
# a single slice with keys = NULL when no slice columns are given.
.slice_index = function(slice_df) {
  n = nrow(slice_df)
  if (ncol(slice_df) == 0) {
    return(list(idx = list(seq_len(n)), keys = NULL))
  }
  keys = interaction(slice_df, drop = TRUE, sep = ".")
  idx = split(seq_len(n), keys)
  keys_df = do.call(rbind, lapply(idx, function(i) slice_df[i[1], , drop = FALSE]))
  list(idx = idx, keys = keys_df)
}

# Moment skewness and excess kurtosis (population moments, documented).
.moments_sk = function(x) {
  x = x[!is.na(x)]
  n = length(x)
  if (n < 2) return(c(skew = NA_real_, kurt = NA_real_))
  m = mean(x)
  s2 = mean((x - m)^2)
  if (s2 <= 0) return(c(skew = NA_real_, kurt = NA_real_))
  c(skew = mean((x - m)^3) / s2^1.5,
    kurt = mean((x - m)^4) / s2^2 - 3)
}

# Attach the stat_infer class and a method label to a result tibble.
.stat_out = function(tbl, label) {
  structure(
    tbl,
    class = c("stat_infer", class(tbl)),
    method_label = label
  )
}

# Two-level factor extraction for group columns (t-test / wilcox / f-test).
.group_two = function(g, n_req = 2L) {
  g = droplevels(as.factor(g))
  if (nlevels(g) != n_req) {
    stop("group must have exactly ", n_req, " levels (got ", nlevels(g),
         ").", call. = FALSE)
  }
  g
}

# Numeric column validation for a selected data frame.
.check_numeric_cols = function(cols_df, arg = ".cols") {
  bad = names(cols_df)[!vapply(cols_df, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(arg, " must select numeric column(s). Non-numeric: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

# Factorize a group column inside a data subset, with level validation.
.to_factor = function(d, group_nm, min_levels = 2L) {
  g = droplevels(as.factor(d[[group_nm]]))
  if (nlevels(g) < min_levels) {
    stop("group must have at least ", min_levels, " levels in the data.",
         call. = FALSE)
  }
  g
}

# ----------------------------------------------------------------------------
# Descriptive statistics
# ----------------------------------------------------------------------------

#' Descriptive Statistics
#'
#' Computes per-variable descriptive statistics (optionally grouped by
#' slices of `.by`): counts, mean, sd, median, IQR, t-based confidence
#' interval for the mean, and moment skewness / excess kurtosis.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric columns.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns;
#'   statistics are computed separately for each slice.
#' @param conf_level Confidence level for the mean CI.
#'
#' @return A tibble with one row per variable (per slice) and class
#'   `stat_infer`: `variable`, `.by` identifiers, `n_used`, `n_missing`,
#'   `mean`, `sd`, `median`, `iqr`, `ci_low`, `ci_high`, `skew`, `kurt`.
#'
#' @examples
#' stat_describe(mtcars, .cols = mpg:hp)
#' stat_describe(mtcars, .cols = mpg, .by = cyl)
#' @export
stat_describe = function(data, .cols, .by = NULL, conf_level = 0.95) {
  .stat_frame(data)
  .check_conf(conf_level)
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  slice_nms = names(slice_df)
  keep = stats::complete.cases(cbind(cols_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  n_drop = nrow(data) - sum(keep)

  sl = .slice_index(slice_df)
  alpha = 1 - conf_level
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cols_df[sl$idx[[k]], , drop = FALSE]
    do.call(rbind, lapply(names(d), function(nm) {
      x = d[[nm]][!is.na(d[[nm]])]
      n = length(x)
      ci_low = ci_high = NA_real_
      if (n >= 2) {
        hw = stats::qt(1 - alpha / 2, df = n - 1) * stats::sd(x) / sqrt(n)
        ci_low = mean(x) - hw
        ci_high = mean(x) + hw
      }
      sk = .moments_sk(x)
      row = tibble::tibble(
        variable = nm,
        n_used = n,
        n_missing = nrow(data) - n,
        mean = if (n) mean(x) else NA_real_,
        sd = if (n >= 2) stats::sd(x) else NA_real_,
        median = if (n) stats::median(x) else NA_real_,
        iqr = if (n) stats::IQR(x) else NA_real_,
        ci_low = ci_low, ci_high = ci_high,
        skew = unname(sk["skew"]), kurt = unname(sk["kurt"])
      )
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(cols_df)
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(slice_nms, setdiff(names(out), slice_nms)), drop = FALSE]
  }
  .stat_out(out, "Descriptive statistics")
}

# ----------------------------------------------------------------------------
# t-test and Wilcoxon tests
# ----------------------------------------------------------------------------

#' t-test (one-sample, two-sample, Welch, paired)
#'
#' Performs t-tests on numeric columns of a tidy data frame. Each test returns
#' one row; multiple columns and `.by` slices are handled by row repetition,
#' with an optional p-value adjustment across rows.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric test
#'   variable(s). One column for one-/two-sample tests; exactly two columns
#'   when `paired = TRUE`.
#' @param group <[`tidy-select`][tidyselect::language]> Optional grouping
#'   column with exactly two levels (the two groups to compare). Required for
#'   two-sample tests, forbidden for paired and one-sample tests.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns;
#'   the test is repeated separately within each slice.
#' @param mu Null value of the mean (or mean difference).
#' @param paired Logical; paired test. Requires exactly two `.cols` and no
#'   `group` (differences are taken row-wise).
#' @param var.equal Logical; pooled-variance t-test (`FALSE` gives Welch).
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`.
#' @param conf_level Confidence level.
#' @param p_adjust P-value adjustment (see [stats::p.adjust]) applied across
#'   the returned rows.
#'
#' @return A tibble with class `stat_infer`: `variable`, `.by` identifiers,
#'   `n_used`, `n_dropped`, `estimate`, `estimate1`, `estimate2`, `statistic`,
#'   `df`, `p.value`, `p.adjusted`, `sig`, `conf_low`, `conf_high`,
#'   `effect` (Cohen's d), `method`, `alternative`.
#'
#' @examples
#' stat_t_test(sleep, .cols = extra, group = group)
#' stat_t_test(mtcars, .cols = mpg, mu = 20)
#' stat_t_test(mtcars, .cols = mpg, group = am, var.equal = TRUE)
#' @export
stat_t_test = function(data, .cols, group, .by = NULL, mu = 0, paired = FALSE,
                       var.equal = FALSE, alternative = "two.sided",
                       conf_level = 0.95, p_adjust = "none") {
  .stat_frame(data)
  alternative = .check_alt(alternative)
  .check_conf(conf_level)
  .check_padjust(p_adjust)
  cols_df = dplyr::select(data, {{ .cols }})
  grp_df = if (missing(group)) {
    dplyr::select(data, dplyr::any_of(character(0)))
  } else {
    dplyr::select(data, {{ group }})
  }
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  overlap = intersect(names(cols_df), c(names(grp_df), names(slice_df)))
  if (length(overlap) > 0) {
    stop(".cols must not overlap group/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  if (ncol(grp_df) > 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  group_nm = if (ncol(grp_df) == 1) names(grp_df)[1] else NULL
  if (paired && (!is.null(group_nm) || ncol(cols_df) != 2)) {
    stop("paired = TRUE requires exactly two .cols and no group.",
         call. = FALSE)
  }
  keep = stats::complete.cases(cbind(cols_df, grp_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  n_drop = nrow(data) - sum(keep)
  has_group = !paired && !is.null(group_nm)
  if (has_group) grp_df[[group_nm]] = .group_two(grp_df[[group_nm]])

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cols_df[sl$idx[[k]], , drop = FALSE]
    g = if (has_group) grp_df[sl$idx[[k]], , drop = FALSE][[group_nm]] else NULL
    vars = if (paired) paste(names(d), collapse = " - ") else names(d)
    do.call(rbind, lapply(seq_along(vars), function(v) {
      fit = if (paired) {
        x = d[[1]]; y = d[[2]]
        stats::t.test(x, y, paired = TRUE, mu = mu,
          alternative = alternative, conf.level = conf_level)
      } else if (has_group) {
        x = d[[v]][g == levels(g)[1]]
        y = d[[v]][g == levels(g)[2]]
        stats::t.test(x, y, mu = mu, var.equal = var.equal,
          alternative = alternative, conf.level = conf_level)
      } else {
        stats::t.test(d[[v]], mu = mu, alternative = alternative,
          conf.level = conf_level)
      }
      effect = if (paired) {
        dv = d[[1]] - d[[2]]
        (mean(dv) - mu) / stats::sd(dv)
      } else if (has_group) {
        x = d[[v]][g == levels(g)[1]]
        y = d[[v]][g == levels(g)[2]]
        n1 = length(x); n2 = length(y)
        sp = sqrt(((n1 - 1) * stats::var(x) + (n2 - 1) * stats::var(y)) /
          (n1 + n2 - 2))
        (mean(x) - mean(y) - mu) / sp
      } else {
        (mean(d[[v]]) - mu) / stats::sd(d[[v]])
      }
      est = fit$estimate
      row = tibble::tibble(
        variable = vars[v],
        n_used = nrow(d),
        n_dropped = n_drop,
        estimate = unname(if (has_group) est[2] - est[1] else est[1]),
        estimate1 = unname(if (has_group) est[1] else NA_real_),
        estimate2 = unname(if (has_group) est[2] else NA_real_),
        statistic = unname(fit$statistic),
        df = unname(fit$parameter),
        p.value = unname(fit$p.value),
        conf_low = unname(fit$conf.int[1]),
        conf_high = unname(fit$conf.int[2]),
        effect = effect,
        method = fit$method,
        alternative = alternative
      )
      if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  out$p.adjusted = stats::p.adjust(out$p.value, method = p_adjust)
  out$sig = .sig_stars(out$p.adjusted)
  .stat_out(out, paste0("t-test (", alternative, ")"))
}

#' Wilcoxon / Mann-Whitney Tests
#'
#' Nonparametric counterparts of the t-test: Wilcoxon signed-rank (one-sample
#' and paired) and Mann-Whitney U (two-sample). Input shapes, `.by` slicing
#' and output layout match [stat_t_test()].
#'
#' @param data A data frame.
#' @inheritParams stat_t_test
#'
#' @return A tibble with class `stat_infer` in the same layout as
#'   [stat_t_test()], with `statistic`, `p.value`, `effect` (rank-biserial
#'   correlation) and no `df` column (NA).
#'
#' @examples
#' stat_wilcox_test(sleep, .cols = extra, group = group)
#' stat_wilcox_test(mtcars, .cols = mpg, mu = 20)
#' @export
stat_wilcox_test = function(data, .cols, group, .by = NULL, mu = 0,
                            paired = FALSE, alternative = "two.sided",
                            conf_level = 0.95, p_adjust = "none") {
  .stat_frame(data)
  alternative = .check_alt(alternative)
  .check_padjust(p_adjust)
  cols_df = dplyr::select(data, {{ .cols }})
  grp_df = if (missing(group)) {
    dplyr::select(data, dplyr::any_of(character(0)))
  } else {
    dplyr::select(data, {{ group }})
  }
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  overlap = intersect(names(cols_df), c(names(grp_df), names(slice_df)))
  if (length(overlap) > 0) {
    stop(".cols must not overlap group/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  if (ncol(grp_df) > 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  group_nm = if (ncol(grp_df) == 1) names(grp_df)[1] else NULL
  if (paired && (!is.null(group_nm) || ncol(cols_df) != 2)) {
    stop("paired = TRUE requires exactly two .cols and no group.",
         call. = FALSE)
  }
  keep = stats::complete.cases(cbind(cols_df, grp_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  n_drop = nrow(data) - sum(keep)
  has_group = !paired && !is.null(group_nm)
  if (has_group) grp_df[[group_nm]] = .group_two(grp_df[[group_nm]])

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cols_df[sl$idx[[k]], , drop = FALSE]
    g = if (has_group) grp_df[sl$idx[[k]], , drop = FALSE][[group_nm]] else NULL
    vars = if (paired) paste(names(d), collapse = " - ") else names(d)
    do.call(rbind, lapply(seq_along(vars), function(v) {
      fit = if (paired) {
        stats::wilcox.test(d[[1]], d[[2]], paired = TRUE, mu = mu,
          alternative = alternative)
      } else if (has_group) {
        stats::wilcox.test(d[[v]] ~ g, mu = mu, alternative = alternative)
      } else {
        stats::wilcox.test(d[[v]], mu = mu, alternative = alternative)
      }
      # Rank-biserial correlation
      effect = if (paired) {
        dv = d[[1]] - d[[2]]
        dv = dv[!is.na(dv)]
        dv = dv[dv != 0]
        n0 = length(dv)
        w_pos = unname(fit$statistic)
        (2 * w_pos - n0 * (n0 + 1) / 2) / (n0 * (n0 + 1) / 2)
      } else if (has_group) {
        n1 = sum(g == levels(g)[1]); n2 = sum(g == levels(g)[2])
        2 * unname(fit$statistic) / (n1 * n2) - 1
      } else {
        NA_real_
      }
      row = tibble::tibble(
        variable = vars[v],
        n_used = nrow(d),
        n_dropped = n_drop,
        estimate = NA_real_,
        estimate1 = NA_real_,
        estimate2 = NA_real_,
        statistic = unname(fit$statistic),
        df = NA_real_,
        p.value = unname(fit$p.value),
        conf_low = NA_real_,
        conf_high = NA_real_,
        effect = effect,
        method = fit$method,
        alternative = alternative
      )
      if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  out$p.adjusted = stats::p.adjust(out$p.value, method = p_adjust)
  out$sig = .sig_stars(out$p.adjusted)
  .stat_out(out, paste0("Wilcoxon rank-sum / signed-rank test (", alternative, ")"))
}

# ----------------------------------------------------------------------------
# ANOVA / Kruskal-Wallis / ANCOVA
# ----------------------------------------------------------------------------

#' Analysis of Variance (ANOVA)
#'
#' Fits an ANOVA on a tidy data frame: `group` columns become factors and are
#' entered with full interactions when several are supplied. One row per
#' model term (plus a `Residuals` row); `resp` and `.by` are repeated
#' (row-per-response and row-per-slice).
#'
#' @param data A data frame.
#' @param resp <[`tidy-select`][tidyselect::language]> Numeric response
#'   column(s); each column is analysed separately.
#' @param group <[`tidy-select`][tidyselect::language]> Factor column(s);
#'   several columns are entered with interactions (`g1 * g2`).
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#'
#' @return A tibble with class `stat_infer`: `resp`, `.by` identifiers,
#'   `term`, `df`, `sumsq`, `meansq`, `statistic`, `p.value`, `p.adjusted`,
#'   `sig`, `effect` (partial eta-squared), `method`.
#'
#' @examples
#' stat_anova(ToothGrowth, resp = len, group = supp)
#' stat_anova(mtcars, resp = mpg, group = c(am, cyl))
#' @export
stat_anova = function(data, resp, group, .by = NULL) {
  .stat_frame(data)
  resp_df = dplyr::select(data, {{ resp }})
  grp_df = dplyr::select(data, {{ group }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(resp_df, arg = "resp")
  if (ncol(grp_df) < 1) {
    stop("group must select at least one column.", call. = FALSE)
  }
  overlap = intersect(names(resp_df), c(names(grp_df), names(slice_df)))
  if (length(overlap) > 0) {
    stop("resp must not overlap group/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(resp_df, grp_df, slice_df))
  resp_df = resp_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cbind(resp_df[sl$idx[[k]], , drop = FALSE],
              grp_df[sl$idx[[k]], , drop = FALSE])
    g_nms = names(grp_df)
    for (nm in g_nms) d[[nm]] = droplevels(as.factor(d[[nm]]))
    if (any(vapply(g_nms, function(nm) nlevels(d[[nm]]) < 2, logical(1)))) {
      stop("each group column needs at least 2 levels in the data.",
           call. = FALSE)
    }
    do.call(rbind, lapply(names(resp_df), function(rnm) {
      rhs = if (length(g_nms) == 1) {
        paste0("`", g_nms, "`")
      } else {
        paste(paste0("`", g_nms, "`"), collapse = " * ")
      }
      f = stats::as.formula(paste0("`", rnm, "` ~ ", rhs))
      fit = stats::aov(f, data = d)
      tab = stats::anova(fit)
      ss_res = tab[nrow(tab), "Sum Sq"]
      term_nms = rownames(tab)
      out = tibble::tibble(
        resp = rnm,
        term = term_nms,
        df = unname(tab$Df),
        sumsq = unname(tab[["Sum Sq"]]),
        meansq = unname(tab[["Mean Sq"]]),
        statistic = unname(tab[["F value"]]),
        p.value = unname(tab[["Pr(>F)"]])
      )
      out$effect = ifelse(term_nms != "Residuals",
        out$sumsq / (out$sumsq + ss_res), NA_real_)
      out
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(resp_df) * length(unique(out$term))
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  out$p.adjusted = stats::p.adjust(out$p.value, method = "none")
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, "Analysis of variance (ANOVA)")
}

#' Kruskal-Wallis Rank Sum Test
#'
#' Nonparametric one-way ANOVA comparing a numeric response across the
#' levels of a single factor column; effect size is epsilon-squared.
#'
#' @param data A data frame.
#' @param resp <[`tidy-select`][tidyselect::language]> Numeric response
#'   column(s); each column is analysed separately.
#' @param group <[`tidy-select`][tidyselect::language]> Single factor column
#'   (two or more levels).
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#'
#' @return A tibble with class `stat_infer`: `resp`, `.by` identifiers,
#'   `n_used`, `statistic`, `df`, `p.value`, `effect` (epsilon-squared),
#'   `method`.
#'
#' @examples
#' stat_kruskal(ToothGrowth, resp = len, group = supp)
#' @export
stat_kruskal = function(data, resp, group, .by = NULL) {
  .stat_frame(data)
  resp_df = dplyr::select(data, {{ resp }})
  grp_df = dplyr::select(data, {{ group }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(resp_df, arg = "resp")
  if (ncol(grp_df) != 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  overlap = intersect(names(resp_df), c(names(grp_df), names(slice_df)))
  if (length(overlap) > 0) {
    stop("resp must not overlap group/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(resp_df, grp_df, slice_df))
  resp_df = resp_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  group_nm = names(grp_df)[1]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cbind(resp_df[sl$idx[[k]], , drop = FALSE],
              grp_df[sl$idx[[k]], , drop = FALSE])
    g = droplevels(as.factor(d[[group_nm]]))
    if (nlevels(g) < 2) {
      stop("group needs at least 2 levels in the data.", call. = FALSE)
    }
    n = nrow(d)
    do.call(rbind, lapply(names(resp_df), function(rnm) {
      fit = stats::kruskal.test(d[[rnm]] ~ g)
      k = nlevels(g)
      tibble::tibble(
        resp = rnm,
        n_used = n,
        statistic = unname(fit$statistic),
        df = unname(fit$parameter),
        p.value = unname(fit$p.value),
        effect = (unname(fit$statistic) - k + 1) / (n - k),
        method = fit$method
      )
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(resp_df)
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, "Kruskal-Wallis rank sum test")
}

#' Analysis of Covariance (ANCOVA)
#'
#' Fits `resp ~ group + covariates` (group as factor, covariates numeric)
#' and returns the ANOVA table with partial eta-squared. Optionally checks
#' the parallel-slopes assumption (group-by-covariate interactions).
#'
#' @param data A data frame.
#' @param resp <[`tidy-select`][tidyselect::language]> Numeric response
#'   column(s); each column is analysed separately.
#' @param group <[`tidy-select`][tidyselect::language]> Single factor column.
#' @param covariates <[`tidy-select`][tidyselect::language]> Numeric
#'   covariate column(s).
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#' @param test_parallel Logical; if `TRUE`, the group-by-covariate
#'   interactions are tested and a warning is issued when any is significant
#'   at 0.05 (parallel-slopes assumption in doubt).
#'
#' @return A tibble with class `stat_infer`: `resp`, `.by` identifiers,
#'   `term`, `df`, `sumsq`, `meansq`, `statistic`, `p.value`, `sig`,
#'   `effect` (partial eta-squared), `method`.
#'
#' @examples
#' d = data.frame(y = c(1:10, 2:11), g = rep(c("a", "b"), each = 10),
#'                x = 1:20)
#' stat_ancova(d, resp = y, group = g, covariates = x)
#' @export
stat_ancova = function(data, resp, group, covariates, .by = NULL,
                       test_parallel = FALSE) {
  .stat_frame(data)
  resp_df = dplyr::select(data, {{ resp }})
  grp_df = dplyr::select(data, {{ group }})
  cov_df = dplyr::select(data, {{ covariates }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(resp_df, arg = "resp")
  .check_numeric_cols(cov_df, arg = "covariates")
  if (ncol(grp_df) != 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  overlap = intersect(names(resp_df), c(names(grp_df), names(cov_df),
                                        names(slice_df)))
  if (length(overlap) > 0) {
    stop("resp must not overlap group/covariates/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(resp_df, grp_df, cov_df, slice_df))
  resp_df = resp_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  cov_df = cov_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  group_nm = names(grp_df)[1]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cbind(resp_df[sl$idx[[k]], , drop = FALSE],
              grp_df[sl$idx[[k]], , drop = FALSE],
              cov_df[sl$idx[[k]], , drop = FALSE])
    d[[group_nm]] = droplevels(as.factor(d[[group_nm]]))
    if (nlevels(d[[group_nm]]) < 2) {
      stop("group needs at least 2 levels in the data.", call. = FALSE)
    }
    cov_nms = names(cov_df)
    if (test_parallel) {
      for (cn in cov_nms) {
        fi = stats::aov(
          stats::as.formula(paste0("`", names(resp_df)[1], "` ~ ",
            paste0("`", group_nm, "` * `", cn, "`"))),
          data = d)
        pi = stats::anova(fi)
        pi = pi[grepl(":", rownames(pi)), , drop = FALSE]
        if (any(pi[["Pr(>F)"]] < 0.05, na.rm = TRUE)) {
          warning("parallel-slopes assumption may not hold: group x `", cn,
                  "` interaction p = ",
                  sprintf("%.4g", min(pi[["Pr(>F)"]], na.rm = TRUE)),
                  ".", call. = FALSE)
        }
      }
    }
    do.call(rbind, lapply(names(resp_df), function(rnm) {
      f = stats::as.formula(paste0("`", rnm, "` ~ `", group_nm, "` + ",
        paste(paste0("`", cov_nms, "`"), collapse = " + ")))
      fit = stats::aov(f, data = d)
      tab = stats::anova(fit)
      ss_res = tab[nrow(tab), "Sum Sq"]
      term_nms = rownames(tab)
      out = tibble::tibble(
        resp = rnm,
        term = term_nms,
        df = unname(tab$Df),
        sumsq = unname(tab[["Sum Sq"]]),
        meansq = unname(tab[["Mean Sq"]]),
        statistic = unname(tab[["F value"]]),
        p.value = unname(tab[["Pr(>F)"]])
      )
      out$effect = ifelse(term_nms != "Residuals",
        out$sumsq / (out$sumsq + ss_res), NA_real_)
      out
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(resp_df) * length(unique(out$term))
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, "Analysis of covariance (ANCOVA)")
}

# ----------------------------------------------------------------------------
# Correlation / chi-square / normality / variance homogeneity
# ----------------------------------------------------------------------------

#' Correlation Analysis
#'
#' Pairwise correlations among numeric columns with test statistics,
#' p-values and (Pearson) confidence intervals.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric columns
#'   (at least 2).
#' @param method `"pearson"`, `"spearman"`, or `"kendall"`. Confidence
#'   intervals are only available for Pearson (`NA` otherwise).
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`.
#' @param conf_level Confidence level (Pearson only).
#' @param p_adjust P-value adjustment across returned rows.
#'
#' @return A tibble with class `stat_infer`: `.by` identifiers, `var1`,
#'   `var2`, `n_used`, `estimate`, `statistic`, `df`, `p.value`,
#'   `p.adjusted`, `sig`, `conf_low`, `conf_high`, `method`.
#'
#' @examples
#' stat_cor(mtcars, .cols = mpg:wt)
#' stat_cor(mtcars, .cols = mpg:wt, method = "spearman")
#' @export
stat_cor = function(data, .cols, method = c("pearson", "spearman", "kendall"),
                    .by = NULL, alternative = "two.sided",
                    conf_level = 0.95, p_adjust = "none") {
  .stat_frame(data)
  method = match.arg(method)
  alternative = .check_alt(alternative)
  .check_conf(conf_level)
  .check_padjust(p_adjust)
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  if (ncol(cols_df) < 2) {
    stop(".cols must select at least 2 columns.", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(cols_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]

  sl = .slice_index(slice_df)
  prs = utils::combn(ncol(cols_df), 2)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cols_df[sl$idx[[k]], , drop = FALSE]
    do.call(rbind, lapply(seq_len(ncol(prs)), function(j) {
      i1 = prs[1, j]; i2 = prs[2, j]
      fit = stats::cor.test(d[[i1]], d[[i2]], method = method,
        alternative = alternative, conf.level = conf_level)
      row = tibble::tibble(
        var1 = names(d)[i1],
        var2 = names(d)[i2],
        n_used = nrow(d),
        estimate = unname(fit$estimate),
        statistic = unname(fit$statistic),
        df = if (method == "pearson") unname(fit$parameter) else NA_real_,
        p.value = unname(fit$p.value),
        conf_low = if (!is.null(fit$conf.int)) unname(fit$conf.int[1]) else NA_real_,
        conf_high = if (!is.null(fit$conf.int)) unname(fit$conf.int[2]) else NA_real_,
        method = fit$method
      )
      if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  out$p.adjusted = stats::p.adjust(out$p.value, method = p_adjust)
  out$sig = .sig_stars(out$p.adjusted)
  .stat_out(out, paste0("Correlation test (", method, ", ", alternative, ")"))
}

#' Chi-square Test of Independence
#'
#' Tests the independence of two factor columns. Reports Cramer's V and
#' warns when expected frequencies are low (suggesting `exact = TRUE`).
#'
#' @param data A data frame.
#' @param row <[`tidy-select`][tidyselect::language]> Single factor column.
#' @param col <[`tidy-select`][tidyselect::language]> Single factor column,
#'   different from `row`.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#' @param correct Logical; continuity correction for 2x2 tables.
#' @param exact Logical; Fisher's exact test instead of chi-square.
#'
#' @return A tibble with class `stat_infer`: `.by` identifiers, `n_used`,
#'   `statistic`, `df`, `p.value`, `sig`, `effect` (Cramer's V),
#'   `min_expected`, `pct_below_5`, `method`.
#'
#' @examples
#' d = data.frame(
#'   treat = rep(c("a", "b"), each = 30),
#'   outcome = sample(c("yes", "no"), 60, TRUE)
#' )
#' stat_chisq(d, row = treat, col = outcome)
#' @export
stat_chisq = function(data, row, col, .by = NULL, correct = TRUE,
                      exact = FALSE) {
  .stat_frame(data)
  row_df = dplyr::select(data, {{ row }})
  col_df = dplyr::select(data, {{ col }})
  slice_df = dplyr::select(data, {{ .by }})
  if (ncol(row_df) != 1) stop("row must select exactly one column.",
    call. = FALSE)
  if (ncol(col_df) != 1) stop("col must select exactly one column.",
    call. = FALSE)
  if (identical(names(row_df), names(col_df))) {
    stop("row and col must be different columns.", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(row_df, col_df, slice_df))
  row_df = row_df[keep, , drop = FALSE]
  col_df = col_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cbind(row_df[sl$idx[[k]], , drop = FALSE],
              col_df[sl$idx[[k]], , drop = FALSE])
    tab = table(as.factor(d[[1]]), as.factor(d[[2]]))
    if (any(dim(tab) < 2)) {
      stop("row and col need at least 2 levels in the data.", call. = FALSE)
    }
    if (exact) {
      fit = stats::fisher.test(tab)
      stat = NA_real_
      df = NA_real_
      method = "Fisher's exact test"
      exp_tab = suppressWarnings(
        stats::chisq.test(tab, correct = FALSE)$expected)
    } else {
      fit = suppressWarnings(stats::chisq.test(tab, correct = correct))
      stat = unname(fit$statistic)
      df = unname(fit$parameter)
      method = fit$method
      exp_tab = fit$expected
    }
    exp_min = min(exp_tab)
    pct5 = mean(exp_tab < 5)
    nr = nrow(tab); nc = ncol(tab)
    chi2 = if (exact) {
      # Cramer's V from the observed chi-square is unavailable for Fisher;
      # compute the uncorrected chi-square statistic for effect size only.
      suppressWarnings(unname(stats::chisq.test(tab, correct = FALSE)$statistic))
    } else {
      suppressWarnings(unname(stats::chisq.test(tab, correct = FALSE)$statistic))
    }
    if (!exact && (pct5 > 0.2 || exp_min < 1)) {
      warning("expected frequencies are low (min = ",
              sprintf("%.2f", exp_min), "; ", sprintf("%.0f%%", 100 * pct5),
              " below 5); consider exact = TRUE.", call. = FALSE)
    }
    row = tibble::tibble(
      n_used = sum(tab),
      statistic = stat,
      df = df,
      p.value = unname(fit$p.value),
      effect = sqrt(chi2 / (sum(tab) * (min(nr, nc) - 1))),
      min_expected = exp_min,
      pct_below_5 = pct5,
      method = method
    )
    if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, "Chi-square test of independence")
}

#' Normality Tests
#'
#' Runs Shapiro-Wilk tests per column (and per `.by` slice), together with
#' moment skewness and excess kurtosis. Shapiro-Wilk requires n between 3
#' and 5000 per group; outside that range the result is `NA`.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric columns.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns
#'   (e.g. groups whose normality is checked separately).
#' @param p_adjust P-value adjustment across returned rows.
#'
#' @return A tibble with class `stat_infer`: `variable`, `.by` identifiers,
#'   `n_used`, `statistic`, `p.value`, `p.adjusted`, `sig`, `skew`, `kurt`,
#'   `method`.
#'
#' @examples
#' stat_normality(mtcars, .cols = mpg:hp)
#' stat_normality(ToothGrowth, .cols = len, .by = supp)
#' @export
stat_normality = function(data, .cols, .by = NULL, p_adjust = "none") {
  .stat_frame(data)
  .check_padjust(p_adjust)
  cols_df = dplyr::select(data, {{ .cols }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  keep = stats::complete.cases(cbind(cols_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cols_df[sl$idx[[k]], , drop = FALSE]
    do.call(rbind, lapply(names(d), function(nm) {
      x = d[[nm]][!is.na(d[[nm]])]
      n = length(x)
      sw = if (n >= 3 && n <= 5000) {
        fit = stats::shapiro.test(x)
        c(stat = unname(fit$statistic), p = unname(fit$p.value))
      } else {
        c(stat = NA_real_, p = NA_real_)
      }
      sk = .moments_sk(x)
      tibble::tibble(
        variable = nm,
        n_used = n,
        statistic = unname(sw["stat"]),
        p.value = unname(sw["p"]),
        skew = unname(sk["skew"]),
        kurt = unname(sk["kurt"]),
        method = "Shapiro-Wilk normality test"
      )
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(cols_df)
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  out$p.adjusted = stats::p.adjust(out$p.value, method = p_adjust)
  out$sig = .sig_stars(out$p.adjusted)
  .stat_out(out, "Shapiro-Wilk normality test")
}

#' Variance Homogeneity Tests
#'
#' Tests the equality of variances across the levels of `group`:
#' Levene (via `car::leveneTest`, centering on median or mean),
#' Bartlett (assumes normality), or the two-level F test.
#'
#' @param data A data frame.
#' @param .cols <[`tidy-select`][tidyselect::language]> Numeric response
#'   column(s); each column is analysed separately.
#' @param group <[`tidy-select`][tidyselect::language]> Single factor column.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice columns.
#' @param method `"levene"`, `"bartlett"`, or `"f"`.
#' @param center For Levene: `"median"` (Brown-Forsythe, default) or
#'   `"mean"`.
#'
#' @return A tibble with class `stat_infer`: `variable`, `.by` identifiers,
#'   `n_used`, `statistic`, `df`, `df2`, `p.value`, `sig`, `method`.
#'
#' @examples
#' stat_hov_test(ToothGrowth, .cols = len, group = supp)
#' stat_hov_test(mtcars, .cols = mpg, group = cyl, method = "bartlett")
#' @export
stat_hov_test = function(data, .cols, group, .by = NULL,
                         method = c("levene", "bartlett", "f"),
                         center = c("median", "mean")) {
  .stat_frame(data)
  method = match.arg(method)
  center = match.arg(center)
  cols_df = dplyr::select(data, {{ .cols }})
  grp_df = dplyr::select(data, {{ group }})
  slice_df = dplyr::select(data, {{ .by }})
  .check_numeric_cols(cols_df)
  if (ncol(grp_df) != 1) {
    stop("group must select exactly one column.", call. = FALSE)
  }
  overlap = intersect(names(cols_df), c(names(grp_df), names(slice_df)))
  if (length(overlap) > 0) {
    stop(".cols must not overlap group/.by columns: ",
         paste(overlap, collapse = ", "), ".", call. = FALSE)
  }
  keep = stats::complete.cases(cbind(cols_df, grp_df, slice_df))
  cols_df = cols_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  group_nm = names(grp_df)[1]

  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = cbind(cols_df[sl$idx[[k]], , drop = FALSE],
              grp_df[sl$idx[[k]], , drop = FALSE])
    g = droplevels(as.factor(d[[group_nm]]))
    if (nlevels(g) < 2) {
      stop("group needs at least 2 levels in the data.", call. = FALSE)
    }
    if (method == "f" && nlevels(g) != 2) {
      stop("method = 'f' requires exactly 2 levels.", call. = FALSE)
    }
    do.call(rbind, lapply(names(cols_df), function(nm) {
      fit = switch(method,
        levene = car::leveneTest(
          stats::reformulate(group_nm, response = nm), data = d,
          center = center),
        bartlett = stats::bartlett.test(
          stats::reformulate(group_nm, response = nm), data = d),
        f = stats::var.test(
          stats::reformulate(group_nm, response = nm), data = d)
      )
      stat = unname(if (method == "levene") fit[1, "F value"] else
        if (method == "bartlett") fit$statistic else fit$statistic)
      df1 = unname(if (method == "levene") fit[1, "Df"] else
        if (method == "bartlett") fit$parameter else fit$parameter[1])
      df2 = if (method == "levene") unname(fit[2, "Df"]) else NA_real_
      p = unname(if (method == "levene") fit[1, "Pr(>F)"] else
        if (method == "bartlett") fit$p.value else fit$p.value)
      mth = switch(method,
        levene = paste0(center, "-centered Levene test"),
        bartlett = "Bartlett test of homogeneity of variances",
        f = "F test to compare two variances")
      tibble::tibble(
        variable = nm,
        n_used = nrow(d),
        statistic = stat,
        df = df1,
        df2 = df2,
        p.value = p,
        method = mth
      )
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_rep = ncol(cols_df)
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_rep), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, "Variance homogeneity test")
}

# ----------------------------------------------------------------------------
# Bootstrap framework (boot-based)
# ----------------------------------------------------------------------------

#' Bootstrap Estimation and Confidence Intervals
#'
#' Generic bootstrap framework built on `boot::boot()`. The `statistic`
#' follows the boot convention `function(data, indices, ...)` — use
#' `indices` to index the original data inside the statistic (the framework
#' never copies the data in the resampling loop). Confidence intervals via
#' `boot::boot.ci()`: BCa (with automatic jackknife), percentile, or basic.
#'
#' @param data A data frame (or any object accepted by `statistic`).
#' @param statistic A function `function(data, indices, ...)` returning a
#'   scalar or a named numeric vector.
#' @param ... Additional arguments passed to `statistic` (and to
#'   `boot::boot()` where unambiguous).
#' @param B Number of bootstrap resamples.
#' @param type `"bca"` (default), `"percentile"`, or `"basic"`.
#' @param conf_level Confidence level.
#' @param seed Optional integer seed for reproducibility.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice
#'   columns; the bootstrap is run separately within each slice.
#' @param parallel `boot::boot` parallel option: `"no"` (default),
#'   `"multicore"` (not on Windows), or `"snow"`.
#' @param ncpus Number of workers for parallel boot.
#'
#' @return A tibble with class `stat_infer`: `.by` identifiers, `term`,
#'   `estimate`, `se`, `bias`, `conf_low`, `conf_high`, `type`, `B`.
#'
#' @examples
#' set.seed(1)
#' stat_bootstrap(mtcars, statistic = \(d, i) mean(d$mpg[i]))
#' stat_bootstrap(mtcars,
#'   statistic = \(d, i) c(mean = mean(d$mpg[i]), sd = stats::sd(d$mpg[i])))
#' @export
stat_bootstrap = function(data, statistic, ..., B = 2000,
                          type = c("bca", "percentile", "basic"),
                          conf_level = 0.95, seed = NULL, .by = NULL,
                          parallel = "no", ncpus = 1) {
  .stat_frame(data)
  type = match.arg(type)
  parallel = match.arg(parallel, c("no", "multicore", "snow"))
  .check_conf(conf_level)
  if (!is.function(statistic)) {
    stop("statistic must be a function of (data, indices, ...).",
         call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1 || B != floor(B)) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  ci_type = switch(type, bca = "bca", percentile = "perc", basic = "basic")
  slice_df = dplyr::select(data, {{ .by }})
  keep = stats::complete.cases(slice_df)
  data = data[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]

  if (!is.null(seed)) set.seed(seed)
  sl = .slice_index(slice_df)
  rows = lapply(seq_along(sl$idx), function(k) {
    d = data[sl$idx[[k]], , drop = FALSE]
    fit = boot::boot(d, statistic, R = B,
      parallel = parallel, ncpus = ncpus, ...)
    t0 = fit$t0
    if (length(t0) == 1) t0 = stats::setNames(t0, "value")
    do.call(rbind, lapply(seq_along(t0), function(j) {
      ci = tryCatch(
        boot::boot.ci(fit, conf = conf_level, type = ci_type, index = j),
        error = function(e) NULL)
      if (is.null(ci)) {
        warning("bootstrap CI (", type, ") failed for component ", j,
                "; returning NA bounds.", call. = FALSE)
        lo = hi = NA_real_
      } else {
        int = ci[[length(ci)]]
        lo = int[4]; hi = int[5]
      }
      t = fit$t[, j]
      tibble::tibble(
        term = if (is.null(names(t0[j])) || is.na(names(t0[j]))) {
          paste0("V", j)
        } else names(t0)[j],
        estimate = unname(t0[j]),
        se = stats::sd(t),
        bias = mean(t) - unname(t0[j]),
        conf_low = lo,
        conf_high = hi,
        type = type,
        B = B
      )
    }))
  })
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    n_comp = length(out$term) / nrow(sl$keys)
    out = dplyr::bind_cols(
      sl$keys[rep(seq_len(nrow(sl$keys)), each = n_comp), , drop = FALSE], out)
    out = out[, c(names(slice_df), setdiff(names(out), names(slice_df))),
              drop = FALSE]
  }
  .stat_out(out, paste0("Bootstrap (", type, " CI, B = ", B, ")"))
}

# ----------------------------------------------------------------------------
# Permutation test framework
# ----------------------------------------------------------------------------

#' Permutation Test
#'
#' Generic permutation-test framework. With `statistic = NULL` the default
#' statistic is chosen by `type`: group mean difference (`two.sample`),
#' mean of differences (`paired`), mean around `mu` (`one.sample`), or
#' Pearson correlation (`independence`). A custom `statistic` of the form
#' `function(data, indices, ...)` is supported for `two.sample` (indices =
#' randomly re-ordered row indices, free permutation) and `independence`
#' (indices = the re-ordered rows used to shuffle the second column).
#'
#' @param data A data frame.
#' @param .cols Numeric test column(s) (required when `statistic = NULL`):
#'   one column for one-/two-sample, two columns for paired/independence.
#' @param group Two-level factor column (required for `two.sample` when
#'   `statistic = NULL`).
#' @param statistic Optional custom function `function(data, indices, ...)`
#'   returning a scalar; observed value uses `indices = seq_len(n)`.
#' @param ... Additional arguments passed to `statistic`.
#' @param type `"two.sample"`, `"paired"`, `"one.sample"`, or
#'   `"independence"`.
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`.
#' @param permutations Number of permutations.
#' @param mu Null value for `one.sample`.
#' @param seed Optional integer seed.
#' @param .by <[`tidy-select`][tidyselect::language]> Optional slice
#'   columns; the test is repeated within each slice.
#'
#' @return A tibble with class `stat_infer`: `.by` identifiers, `observed`,
#'   `p.value`, `null_mean`, `null_sd`, `permutations`, `type`, `alternative`.
#'
#' @examples
#' set.seed(1)
#' d = data.frame(y = c(rnorm(20), rnorm(20, 0.8)),
#'                g = rep(c("a", "b"), each = 20))
#' stat_perm_test(d, .cols = y, group = g, type = "two.sample")
#' @export
stat_perm_test = function(data, .cols = NULL, group = NULL, statistic = NULL,
                          ..., type = c("two.sample", "paired", "one.sample",
                                       "independence"),
                          alternative = "two.sided", permutations = 2000,
                          mu = 0, seed = NULL, .by = NULL) {
  .stat_frame(data)
  type = match.arg(type)
  alternative = .check_alt(alternative)
  if (!is.numeric(permutations) || length(permutations) != 1 ||
      permutations < 1 || permutations != floor(permutations)) {
    stop("permutations must be a positive integer.", call. = FALSE)
  }
  cols_df = if (missing(.cols)) {
    dplyr::select(data, dplyr::any_of(character(0)))
  } else {
    dplyr::select(data, {{ .cols }})
  }
  grp_df = if (missing(group)) {
    dplyr::select(data, dplyr::any_of(character(0)))
  } else {
    dplyr::select(data, {{ group }})
  }
  slice_df = dplyr::select(data, {{ .by }})
  slice_nms = names(slice_df)

  .run_one = function(d, cols_df, grp_df) {
    B = permutations
    if (!is.null(statistic)) {
      if (!is.function(statistic)) {
        stop("statistic must be a function of (data, indices, ...).",
             call. = FALSE)
      }
      if (!type %in% c("two.sample", "independence")) {
        stop("custom statistic supports type 'two.sample' and ",
             "'independence' only.", call. = FALSE)
      }
      n = nrow(d)
      T0 = statistic(d, seq_len(n), ...)
      if (type == "two.sample") {
        Tstar = vapply(seq_len(B), function(b) {
          statistic(d, sample.int(n), ...)
        }, numeric(1))
      } else {
        idx_mat = replicate(B, sample.int(n))
        Tstar = vapply(seq_len(B), function(b) {
          d2 = d
          d2[[ncol(d)]] = d[[ncol(d)]][idx_mat[, b]]
          statistic(d2, seq_len(n), ...)
        }, numeric(1))
      }
    } else {
      .check_numeric_cols(cols_df)
      if (type == "two.sample") {
        if (ncol(cols_df) != 1 || ncol(grp_df) != 1) {
          stop("two.sample requires one .cols column and a group column.",
               call. = FALSE)
        }
        g = droplevels(as.factor(grp_df[[1]]))
        if (nlevels(g) != 2) {
          stop("group must have exactly 2 levels.", call. = FALSE)
        }
        x = cols_df[[1]][g == levels(g)[1]]
        y = cols_df[[1]][g == levels(g)[2]]
        T0 = mean(x) - mean(y)
        lab = rep(c(FALSE, TRUE), c(length(x), length(y)))
        Tstar = vapply(seq_len(B), function(b) {
          lp = sample(lab)
          mean(c(x, y)[lp]) - mean(c(x, y)[!lp])
        }, numeric(1))
      } else if (type == "paired") {
        if (ncol(cols_df) != 2) {
          stop("paired requires exactly two .cols columns.", call. = FALSE)
        }
        dv = cols_df[[1]] - cols_df[[2]]
        T0 = mean(dv)
        sgn = matrix(sample(c(-1, 1), length(dv) * B, replace = TRUE),
          nrow = B)
        Tstar = as.vector(sgn %*% dv) / length(dv)
      } else if (type == "one.sample") {
        if (ncol(cols_df) != 1) {
          stop("one.sample requires exactly one .cols column.", call. = FALSE)
        }
        dv = cols_df[[1]] - mu
        T0 = mean(dv)
        sgn = matrix(sample(c(-1, 1), length(dv) * B, replace = TRUE),
          nrow = B)
        Tstar = as.vector(sgn %*% dv) / length(dv)
      } else {
        if (ncol(cols_df) != 2) {
          stop("independence requires exactly two .cols columns.",
               call. = FALSE)
        }
        x = cols_df[[1]]; y = cols_df[[2]]
        T0 = stats::cor(x, y)
        idx_mat = replicate(B, sample.int(length(y)))
        Tstar = apply(idx_mat, 2, function(i) stats::cor(x, y[i]))
      }
    }
    if (alternative == "greater") {
      p = (1 + sum(Tstar >= T0)) / (B + 1)
    } else if (alternative == "less") {
      p = (1 + sum(Tstar <= T0)) / (B + 1)
    } else {
      c0 = if (!is.null(statistic)) mean(Tstar) else 0
      p = (1 + sum(abs(Tstar - c0) >= abs(T0 - c0))) / (B + 1)
    }
    tibble::tibble(
      observed = unname(T0),
      p.value = p,
      null_mean = mean(Tstar),
      null_sd = stats::sd(Tstar),
      permutations = B,
      type = type,
      alternative = alternative
    )
  }

  # Complete rows across the columns actually used
  used_cols = names(cols_df)
  keep = stats::complete.cases(cbind(cols_df, grp_df, slice_df))
  data = data[keep, , drop = FALSE]
  cols_df = cols_df[keep, , drop = FALSE]
  grp_df = grp_df[keep, , drop = FALSE]
  slice_df = slice_df[keep, , drop = FALSE]
  # If a custom statistic is used, operate on the full (kept) data frame
  if (!is.null(statistic)) {
    sl = .slice_index(slice_df)
    rows = lapply(seq_along(sl$idx), function(k) {
      d = data[sl$idx[[k]], , drop = FALSE]
      row = .run_one(d, NULL, NULL)
      if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
    })
  } else {
    sl = .slice_index(slice_df)
    rows = lapply(seq_along(sl$idx), function(k) {
      i = sl$idx[[k]]
      row = .run_one(data[i, , drop = FALSE],
        cols_df[i, , drop = FALSE], grp_df[i, , drop = FALSE])
      if (is.null(sl$keys)) row else dplyr::bind_cols(sl$keys[k, , drop = FALSE], row)
    })
  }
  out = tibble::as_tibble(do.call(rbind, rows))
  if (!is.null(sl$keys)) {
    out = out[, c(slice_nms, setdiff(names(out), slice_nms)), drop = FALSE]
  }
  out$sig = .sig_stars(out$p.value)
  .stat_out(out, paste0("Permutation test (", type, ", ", alternative, ")"))
}

# ----------------------------------------------------------------------------
# Print method
# ----------------------------------------------------------------------------

#' @export
print.stat_infer = function(x, ...) {
  cat(attr(x, "method_label"), "\n")
  class(x) = setdiff(class(x), "stat_infer")
  NextMethod()
}
