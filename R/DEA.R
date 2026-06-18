# ============================================================================
# DEA Models: Pure R implementation using lpSolveAPI
# ============================================================================
# Supports: CCR, BCC, SBM, super-efficiency, undesirable outputs, Malmquist
# Reference: Charnes, Cooper & Rhodes (1978); Banker, Charnes & Cooper (1984);
#            Tone (2001); Andersen & Petersen (1993); Färe et al. (1994)
# ============================================================================

# ============================================================================
# Internal helpers: Data extraction and validation
# ============================================================================

#' Extract DMU names and matrices from a data frame
#' Returns matrices in m×n format (inputs/outputs as rows, DMUs as columns),
#' matching the MATLAB convention used in DEA linear programming.
#' @noRd
extract_dea_data = function(data, inputs, outputs, ud_outputs = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (ncol(data) < 3) stop("`data` must have at least 3 columns (DMU + >=1 input + >=1 output).")

  dmu_names = data[[1]]
  X = t(as.matrix(data[, inputs, drop = FALSE]))      # m × n
  Y_all = t(as.matrix(data[, outputs, drop = FALSE]))   # s × n

  if (anyNA(X) || anyNA(Y_all)) stop("Input/output columns must not contain NA values.")

  if (!is.null(ud_outputs)) {
    ud_idx = ud_outputs
    if (max(ud_idx) > nrow(Y_all) || min(ud_idx) < 1)
      stop("`ud_outputs` must be valid indices within `outputs`.")
    Y_g = Y_all[-ud_idx, , drop = FALSE]   # desirable outputs
    Y_b = Y_all[ ud_idx, , drop = FALSE]   # undesirable outputs
  } else {
    Y_g = Y_all
    Y_b = NULL
  }

  if (nrow(Y_g) == 0) stop("At least one desirable output is required.")

  list(dmu = dmu_names,
       X = X, Y_g = Y_g, Y_b = Y_b,
       n = ncol(X), m = nrow(X), s1 = nrow(Y_g),
       s2 = if (is.null(Y_b)) 0L else nrow(Y_b))
}

# ============================================================================
# Internal helpers: LP solvers
# ============================================================================

#' Solve radial DEA LP (CCR/BCC) for a single DMU
#'
#' Matrices MUST be in m×n format (inputs/outputs as rows, DMUs as columns).
#'
#' @param X Input matrix (m × n_ref)
#' @param Y Desirable output matrix (s1 × n_ref)
#' @param x0 Target DMU input vector (length m)
#' @param y0 Target DMU desirable output vector (length s1)
#' @param orientation "io" or "oo"
#' @param rts "crs" or "vrs"
#' @return list(theta, phi, lambda, s_in, s_out)
#'   theta = Farrell efficiency (0~1); phi = Shephard output distance (≥1, NA for io)
#' @noRd
solve_radial_lp = function(X, Y, x0, y0, orientation = "io", rts = "vrs") {
  n_ref = ncol(X)
  m = length(x0)
  q = length(y0)

  lp = lpSolveAPI::make.lp(0, n_ref + 1)

  if (orientation == "io") {
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), 1))  # min θ
    for (k in seq_len(m))
      lpSolveAPI::add.constraint(lp, c(X[k, ], -x0[k]), "<=", 0)
    for (k in seq_len(q))
      lpSolveAPI::add.constraint(lp, c(-Y[k, ], 0), "<=", -y0[k])
  } else {
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), -1))  # min -φ ≡ max φ
    for (k in seq_len(m))
      lpSolveAPI::add.constraint(lp, c(X[k, ], 0), "<=", x0[k])
    for (k in seq_len(q))
      lpSolveAPI::add.constraint(lp, c(-Y[k, ], y0[k]), "<=", 0)
  }

  if (rts == "vrs")
    lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "=", 1)

  lpSolveAPI::set.bounds(lp, lower = rep(0, n_ref + 1))

  status = lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(theta = NA_real_, phi = NA_real_,
                lambda = rep(NA_real_, n_ref),
                s_in = rep(NA_real_, m), s_out = rep(NA_real_, q)))
  }

  sol = lpSolveAPI::get.variables(lp)
  lambda = sol[1:n_ref]
  theta_raw = sol[n_ref + 1]

  if (orientation == "io") {
    s_in  = pmax(0, theta_raw * x0 - as.vector(X %*% lambda))
    s_out = pmax(0, as.vector(Y %*% lambda) - y0)
    eff = theta_raw
    phi = NA_real_
  } else {
    phi = theta_raw
    s_in  = pmax(0, x0 - as.vector(X %*% lambda))
    s_out = pmax(0, as.vector(Y %*% lambda) - phi * y0)
    eff = 1 / phi
  }

  list(theta = eff, phi = phi,
       lambda = lambda, s_in = s_in, s_out = s_out)
}

#' Solve SBM LP (Tone 2001, Charnes-Cooper transformed) for a single DMU
#'
#' @param X Input matrix (m × n_ref)
#' @param Y_g Desirable output matrix (s1 × n_ref)
#' @param Y_b Undesirable output matrix (s2 × n_ref), NULL if none
#' @param x0 Target DMU input vector (length m)
#' @param y0_g Target DMU desirable output vector (length s1)
#' @param y0_b Target DMU undesirable output vector (length s2), NULL if none
#' @param rts "crs" or "vrs"
#' @return list(rho, s_in, s_g, s_b, lambda)
#' @noRd
solve_sbm_lp = function(X, Y_g, Y_b, x0, y0_g, y0_b, rts = "vrs") {
  n_ref = ncol(X)
  m = length(x0)
  s1 = length(y0_g)
  s2 = if (is.null(Y_b)) 0L else nrow(Y_b)

  nvar = m + s1 + s2 + n_ref + 1
  c_val = 1 / max(s1 + s2, 1)

  lp = lpSolveAPI::make.lp(0, nvar)

  obj = numeric(nvar)
  obj[1:m] = -1 / (m * x0)
  obj[nvar] = 1
  lpSolveAPI::set.objfn(lp, obj)

  # Constraint 1: output slack weighted sum + t = 1
  con1 = numeric(nvar)
  con1[(m + 1):(m + s1)] = c_val / y0_g
  if (s2 > 0) con1[(m + s1 + 1):(m + s1 + s2)] = c_val / y0_b
  con1[nvar] = 1
  lpSolveAPI::add.constraint(lp, con1, "=", 1)

  # Input slack constraints: s^- + X λ - t x0 = 0
  ii_lambda = (m + s1 + s2 + 1):(m + s1 + s2 + n_ref)
  for (k in seq_len(m)) {
    row = numeric(nvar)
    row[k] = 1
    row[ii_lambda] = X[k, ]
    row[nvar] = -x0[k]
    lpSolveAPI::add.constraint(lp, row, "=", 0)
  }

  # Desirable output slack constraints: -s^g + Y_g λ - t y0_g = 0
  for (k in seq_len(s1)) {
    row = numeric(nvar)
    row[m + k] = -1
    row[ii_lambda] = Y_g[k, ]
    row[nvar] = -y0_g[k]
    lpSolveAPI::add.constraint(lp, row, "=", 0)
  }

  # Undesirable output slack constraints: s^b + Y_b λ - t y0_b = 0
  if (s2 > 0) {
    for (k in seq_len(s2)) {
      row = numeric(nvar)
      row[m + s1 + k] = 1
      row[ii_lambda] = Y_b[k, ]
      row[nvar] = -y0_b[k]
      lpSolveAPI::add.constraint(lp, row, "=", 0)
    }
  }

  # VRS constraint: Σλ = t
  if (rts == "vrs") {
    row = numeric(nvar)
    row[ii_lambda] = 1
    row[nvar] = -1
    lpSolveAPI::add.constraint(lp, row, "=", 0)
  }

  lpSolveAPI::set.bounds(lp, lower = rep(0, nvar))

  status = lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(rho = NA_real_, s_in = rep(NA_real_, m),
                s_g = rep(NA_real_, s1), s_b = rep(NA_real_, max(s2, 0)),
                lambda = rep(NA_real_, n_ref)))
  }

  sol = lpSolveAPI::get.variables(lp)
  t_val = sol[nvar]

  rho = lpSolveAPI::get.objective(lp)

  s_in  = sol[1:m] / t_val
  s_g   = sol[(m + 1):(m + s1)] / t_val
  s_b   = if (s2 > 0) sol[(m + s1 + 1):(m + s1 + s2)] / t_val else numeric(0)
  lambda = sol[ii_lambda] / t_val

  list(rho = rho, s_in = s_in, s_g = s_g, s_b = s_b, lambda = lambda)
}

# ============================================================================
# Internal helpers: Result assembly
# ============================================================================

#' Build the output list for DEA functions
#' @noRd
build_result = function(dd, theta_vec, lambda_mat, s_in_mat, s_out_mat, rts) {
  rts_label = toupper(rts)

  list(
    efficiencies = theta_vec,
    lambdas = lambda_mat,
    slacks = list(
      inputs   = s_in_mat,
      outputs  = s_out_mat
    ),
    targets = list(
      inputs   = dd$X - t(s_in_mat),
      outputs  = dd$Y_g + t(s_out_mat[, seq_len(dd$s1), drop = FALSE])
    ),
    returns = structure(rep(rts_label, dd$n), names = dd$dmu),
    model = "DEA",
    orientation = attr(dd, "orientation") %||% "io",
    dmu = dd$dmu
  )
}

# ============================================================================
# Public DEA functions
# ============================================================================

#' DEA efficiency analysis
#'
#' @name DEA
#' @param data A data frame. 1st column = DMU names.
#' @param inputs Column indices/names of input variables.
#' @param outputs Column indices/names of output variables.
#' @param ud_outputs Column indices (within `outputs`) of undesirable outputs, or NULL.
#' @param orientation `"io"` (input-oriented, default) or `"oo"` (output-oriented).
#'   All DEA functions return Farrell efficiency (0~1), where 1 = efficient.
#'   For output orientation: value = 1/φ, where φ ≥ 1 is the Shephard distance.
#' @param rts `"vrs"` (variable returns to scale, default) or `"crs"` (constant).
#' @param x0,y0 Used by Malmquist cross-period evaluation (not for direct calls).
#' @return A list with components: `efficiencies`, `lambdas`, `slacks`,
#'   `targets`, `returns`, `model`, `orientation`, `dmu`.
NULL

#' @describeIn DEA Standard radial DEA (CCR/BCC)
#' @examples
#' df = data.frame(
#'   DMU = paste0("DMU", 1:7),
#'   x1  = c(20, 60, 40, 60, 70, 30, 50),
#'   x2  = c(151, 200, 120, 170, 250, 210, 90),
#'   y1  = c(100, 210, 150, 240, 220, 80, 200)
#' )
#' basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
#' basic_DEA(df, inputs = 2:3, outputs = 4, rts = "vrs")
#' @export
basic_DEA = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  orientation = match.arg(orientation, c("io", "oo"))
  rts = match.arg(rts, c("vrs", "crs"))

  dd = extract_dea_data(data, inputs, outputs, ud_outputs)
  n = dd$n; m = dd$m; s1 = dd$s1; s2 = dd$s2

  theta_vec = setNames(rep(NA_real_, n), dd$dmu)
  lambda_mat = matrix(NA_real_, n, n,
    dimnames = list(dd$dmu, dd$dmu))
  s_in_mat  = matrix(NA_real_, n, m,
    dimnames = list(dd$dmu, colnames(t(dd$X))))
  s_out_mat = matrix(NA_real_, n, s1 + s2,
    dimnames = list(dd$dmu,
      c(colnames(t(dd$Y_g)),
        if (s2 > 0) colnames(t(dd$Y_b)) else NULL)))

  for (i in seq_len(n)) {
    y_g = dd$Y_g[, i]
    y_b = if (s2 > 0) dd$Y_b[, i] else NULL
    res = solve_radial_lp(X = dd$X, Y = dd$Y_g,
                          x0 = dd$X[, i], y0 = y_g,
                          orientation = orientation, rts = rts)
    if (!is.na(res$theta)) {
      theta_vec[i] = res$theta
      lambda_mat[i, ] = res$lambda
      s_in_mat[i, ] = res$s_in
      s_out_mat[i, seq_len(s1)] = res$s_out
      if (s2 > 0) s_out_mat[i, (s1 + 1):(s1 + s2)] = rep(NA_real_, s2)
    }
  }

  attr(dd, "orientation") = orientation
  build_result(dd, theta_vec, lambda_mat, s_in_mat, s_out_mat, rts)
}

#' @describeIn DEA Super-efficiency DEA (self excluded, radial only)
#' @examples
#' df = data.frame(
#'   DMU = paste0("DMU", 1:7),
#'   x1  = c(20, 60, 40, 60, 70, 30, 50),
#'   x2  = c(151, 200, 120, 170, 250, 210, 90),
#'   y1  = c(100, 210, 150, 240, 220, 80, 200)
#' )
#' super_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
#' @export
super_DEA = function(data, inputs, outputs,
                     orientation = "io", rts = "vrs") {
  orientation = match.arg(orientation, c("io", "oo"))
  rts = match.arg(rts, c("vrs", "crs"))

  dd = extract_dea_data(data, inputs, outputs, ud_outputs = NULL)
  n = dd$n; m = dd$m; s1 = dd$s1

  theta_vec = setNames(rep(NA_real_, n), dd$dmu)
  lambda_mat = matrix(NA_real_, n, n,
    dimnames = list(dd$dmu, dd$dmu))
  s_in_mat  = matrix(NA_real_, n, m,
    dimnames = list(dd$dmu, colnames(t(dd$X))))
  s_out_mat = matrix(NA_real_, n, s1,
    dimnames = list(dd$dmu, colnames(t(dd$Y_g))))

  for (i in seq_len(n)) {
    if (n <= 1) next
    idx_ref = setdiff(seq_len(n), i)
    res = solve_radial_lp(
      X = dd$X[, idx_ref, drop = FALSE],
      Y = dd$Y_g[, idx_ref, drop = FALSE],
      x0 = dd$X[, i], y0 = dd$Y_g[, i],
      orientation = orientation, rts = rts)
    if (!is.na(res$theta)) {
      theta_vec[i] = res$theta
      lambda_mat[i, idx_ref] = res$lambda
      s_in_mat[i, ]  = res$s_in
      s_out_mat[i, seq_len(s1)] = res$s_out
    }
  }

  attr(dd, "orientation") = orientation
  build_result(dd, theta_vec, lambda_mat, s_in_mat, s_out_mat, rts)
}

#' @describeIn DEA Standard Slacks-Based Measure (SBM, Tone 2001)
#' @examples
#' df = data.frame(
#'   DMU = paste0("DMU", 1:7),
#'   x1  = c(20, 60, 40, 60, 70, 30, 50),
#'   x2  = c(151, 200, 120, 170, 250, 210, 90),
#'   y1  = c(100, 210, 150, 240, 220, 80, 200)
#' )
#' basic_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
#' @export
basic_SBM = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  rts = match.arg(rts, c("vrs", "crs"))

  dd = extract_dea_data(data, inputs, outputs, ud_outputs)
  n = dd$n; m = dd$m; s1 = dd$s1; s2 = dd$s2

  rho_vec = setNames(rep(NA_real_, n), dd$dmu)
  lambda_mat = matrix(NA_real_, n, n,
    dimnames = list(dd$dmu, dd$dmu))
  s_in_mat  = matrix(NA_real_, n, m,
    dimnames = list(dd$dmu, colnames(t(dd$X))))
  s_out_mat = matrix(NA_real_, n, s1 + s2,
    dimnames = list(dd$dmu,
      c(colnames(t(dd$Y_g)),
        if (s2 > 0) colnames(t(dd$Y_b)) else NULL)))

  for (i in seq_len(n)) {
    res = solve_sbm_lp(
      X = dd$X, Y_g = dd$Y_g, Y_b = dd$Y_b,
      x0 = dd$X[, i], y0_g = dd$Y_g[, i],
      y0_b = if (s2 > 0) dd$Y_b[, i] else NULL,
      rts = rts)
    if (!is.na(res$rho)) {
      rho_vec[i] = res$rho
      lambda_mat[i, ] = res$lambda
      s_in_mat[i, ] = res$s_in
      s_out_mat[i, seq_len(s1)] = res$s_g
      if (s2 > 0) s_out_mat[i, (s1 + 1):(s1 + s2)] = res$s_b
    }
  }

  build_result(dd, rho_vec, lambda_mat, s_in_mat, s_out_mat, rts)
}

#' @describeIn DEA Super-efficiency SBM (self excluded, no undesirable outputs)
#' @examples
#' df = data.frame(
#'   DMU = paste0("DMU", 1:7),
#'   x1  = c(20, 60, 40, 60, 70, 30, 50),
#'   x2  = c(151, 200, 120, 170, 250, 210, 90),
#'   y1  = c(100, 210, 150, 240, 220, 80, 200)
#' )
#' super_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
#' @export
super_SBM = function(data, inputs, outputs,
                     orientation = "io", rts = "vrs") {
  rts = match.arg(rts, c("vrs", "crs"))

  dd = extract_dea_data(data, inputs, outputs, ud_outputs = NULL)
  n = dd$n; m = dd$m; s1 = dd$s1

  rho_vec = setNames(rep(NA_real_, n), dd$dmu)
  lambda_mat = matrix(NA_real_, n, n,
    dimnames = list(dd$dmu, dd$dmu))
  s_in_mat  = matrix(NA_real_, n, m,
    dimnames = list(dd$dmu, colnames(t(dd$X))))
  s_out_mat = matrix(NA_real_, n, s1,
    dimnames = list(dd$dmu, colnames(t(dd$Y_g))))

  for (i in seq_len(n)) {
    if (n <= 1) next
    idx_ref = setdiff(seq_len(n), i)
    res = solve_sbm_lp(
      X = dd$X[, idx_ref, drop = FALSE],
      Y_g = dd$Y_g[, idx_ref, drop = FALSE],
      Y_b = NULL,
      x0 = dd$X[, i], y0_g = dd$Y_g[, i], y0_b = NULL,
      rts = rts)
    if (!is.na(res$rho)) {
      rho_vec[i] = res$rho
      lambda_mat[i, idx_ref] = res$lambda
      s_in_mat[i, ]  = res$s_in
      s_out_mat[i, ] = res$s_g
    }
  }

  build_result(dd, rho_vec, lambda_mat, s_in_mat, s_out_mat, rts)
}

# ============================================================================
# Malmquist Productivity Index
# ============================================================================


#' @describeIn DEA Malmquist productivity index
#' @examples
#' panel = data.frame(
#'   DMU    = rep(paste0("DMU", 1:5), 3),
#'   Period = rep(1:3, each = 5),
#'   x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
#'   y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
#' )
#' malmquist(panel, period = "Period", inputs = 3, outputs = 4,
#'   rts = "crs", type1 = "cont", type2 = "fgnz")
#' @export
malmquist = function(data, period, inputs, outputs,
                     orientation = "oo", rts = "vrs",
                     type1 = "glob", type2 = "rd") {
  orientation = match.arg(orientation, c("io", "oo"))
  rts = match.arg(rts, c("vrs", "crs"))
  type1 = match.arg(type1, c("cont", "seq", "glob"))
  type2 = match.arg(type2, c("fgnz", "rd"))

  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (ncol(data) < 4)       stop("data must have >= 4 columns.")

  dmu_vec    = data[[1]]
  period_vec = data[[period]]
  X_all = t(as.matrix(data[, inputs,  drop = FALSE]))
  Y_all = t(as.matrix(data[, outputs, drop = FALSE]))

  periods = sort(unique(period_vec))
  nt = length(periods)
  if (nt < 2) stop("At least 2 periods required.")

  dmu0 = dmu_vec[period_vec == periods[1]]
  nd   = length(dmu0)

  # ---- 1. Contemporaneous Shephard distances ----
  D_crs = matrix(NA_real_, nt, nd, dimnames = list(periods, dmu0))
  D_vrs = matrix(NA_real_, nt, nd, dimnames = list(periods, dmu0))

  for (t_idx in seq_len(nt)) {
    idx = which(period_vec == periods[t_idx])
    dmu_t = dmu_vec[idx]
    X_ref = X_all[, idx, drop = FALSE]
    Y_ref = Y_all[, idx, drop = FALSE]
    for (j in seq_len(nd)) {
      ii = match(dmu0[j], dmu_t)
      if (is.na(ii)) next
      x0 = X_all[, idx[ii]]; y0 = Y_all[, idx[ii]]
      D_crs[t_idx, j] = dist_shephard(x0, y0, X_ref, Y_ref, orientation, "crs")
      D_vrs[t_idx, j] = dist_shephard(x0, y0, X_ref, Y_ref, orientation, "vrs")
    }
  }

  # ---- 2. Global Shephard distances ----
  Dg_crs = matrix(NA_real_, nt, nd, dimnames = list(periods, dmu0))
  Dg_vrs = matrix(NA_real_, nt, nd, dimnames = list(periods, dmu0))

  if (type1 == "glob") {
    for (t_idx in seq_len(nt)) {
      idx = which(period_vec == periods[t_idx])
      dmu_t = dmu_vec[idx]
      for (j in seq_len(nd)) {
        ii = match(dmu0[j], dmu_t)
        if (is.na(ii)) next
        x0 = X_all[, idx[ii]]; y0 = Y_all[, idx[ii]]
        Dg_crs[t_idx, j] = dist_shephard(x0, y0, X_all, Y_all, orientation, "crs")
        Dg_vrs[t_idx, j] = dist_shephard(x0, y0, X_all, Y_all, orientation, "vrs")
      }
    }
  }

  # ---- 3. Cross-period Shephard distances ----
  D12_crs = matrix(NA_real_, nt - 1, nd, dimnames = list(periods[-1], dmu0))
  D12_vrs = matrix(NA_real_, nt - 1, nd, dimnames = list(periods[-1], dmu0))
  D21_crs = matrix(NA_real_, nt - 1, nd, dimnames = list(periods[-1], dmu0))

  for (t_idx in seq_len(nt - 1)) {
    p_t   = periods[t_idx]
    p_t1  = periods[t_idx + 1]
    idx_t  = which(period_vec == p_t)
    idx_t1 = which(period_vec == p_t1)
    X_t  = X_all[, idx_t,  drop = FALSE]
    Y_t  = Y_all[, idx_t,  drop = FALSE]
    X_t1 = X_all[, idx_t1, drop = FALSE]
    Y_t1 = Y_all[, idx_t1, drop = FALSE]

    for (j in seq_len(nd)) {
      ii_t  = match(dmu0[j], dmu_vec[idx_t])
      ii_t1 = match(dmu0[j], dmu_vec[idx_t1])
      if (is.na(ii_t) || is.na(ii_t1)) next

      x_t  = X_all[, idx_t[ii_t]];   y_t  = Y_all[, idx_t[ii_t]]
      x_t1 = X_all[, idx_t1[ii_t1]]; y_t1 = Y_all[, idx_t1[ii_t1]]

      D12_crs[t_idx, j] = dist_shephard(x_t1, y_t1, X_t,  Y_t,  orientation, "crs")
      D12_vrs[t_idx, j] = dist_shephard(x_t1, y_t1, X_t,  Y_t,  orientation, "vrs")
      D21_crs[t_idx, j] = dist_shephard(x_t,  y_t,  X_t1, Y_t1, orientation, "crs")
    }
  }

  # ---- 4. Convert Shephard → Farrell ----
  # deaR uses Farrell efficiencies internally to compute indices.
  # Input: Farrell = θ. Output: Farrell = 1/φ.
  to_farrell = if (orientation == "oo") function(d) 1 / d else identity

  E    = to_farrell(D_crs);    Ev   = to_farrell(D_vrs)
  Eg   = to_farrell(Dg_crs);   Egv  = to_farrell(Dg_vrs)
  E12  = to_farrell(D12_crs);  Ev12 = to_farrell(D12_vrs)
  E21  = to_farrell(D21_crs)

  # ---- 5. Compute indices (deaR formulas) ----
  out = vector("list", (nt - 1) * nd)
  k  = 0L

  for (t_idx in seq_len(nt - 1)) {
    label = sprintf("%s~%s", periods[t_idx], periods[t_idx + 1])
    for (j in seq_len(nd)) {

      if (type1 == "glob") {
        if (rts == "crs") {
          ec_i = E[t_idx + 1, j] / E[t_idx, j]
          tc_i = (Eg[t_idx + 1, j] * E[t_idx, j]) /
                 (E[t_idx + 1, j] * Eg[t_idx, j])
          mi_i = ec_i * tc_i
          pech_i = NA_real_; sech_i = NA_real_
        } else {
          pech_i = Ev[t_idx + 1, j] / Ev[t_idx, j]
          tc_i   = (Egv[t_idx + 1, j] * Ev[t_idx, j]) /
                   (Ev[t_idx + 1, j] * Egv[t_idx, j])
          if (type2 == "rd") {
            se_i_t     = Ev[t_idx, j]  / E[t_idx, j]
            se_i_cross = Ev12[t_idx, j] / E12[t_idx, j]
            sech_i = se_i_t / se_i_cross
          } else {
            sech_i = NA_real_
          }
          mi_i = tc_i * pech_i * sech_i
          ec_i = pech_i * sech_i
        }

      } else {
        # cont / seq: Färe et al. (1994) geometric mean
        tc_i = sqrt((E12[t_idx, j] / E[t_idx + 1, j]) *
                    (E[t_idx, j]    / E21[t_idx, j]))
        if (rts == "crs") {
          ec_i = E[t_idx + 1, j] / E[t_idx, j]
          mi_i = ec_i * tc_i
          pech_i = NA_real_; sech_i = NA_real_
        } else {
          pech_i = Ev[t_idx + 1, j] / Ev[t_idx, j]
          if (type2 == "rd") {
            se_i_t     = Ev[t_idx, j]  / E[t_idx, j]
            se_i_cross = Ev12[t_idx, j] / E12[t_idx, j]
            sech_i = se_i_t / se_i_cross
          } else {
            se_i_t  = Ev[t_idx, j]  / E[t_idx, j]
            se_i_t1 = Ev[t_idx + 1, j] / E[t_idx + 1, j]
            sech_i = se_i_t / se_i_t1
          }
          mi_i = tc_i * pech_i * sech_i
          ec_i = pech_i * sech_i
        }
      }

      k = k + 1L
      out[[k]] = data.frame(
        Period = label, DMU = dmu0[j],
        mi = mi_i, ec = ec_i, tc = tc_i,
        pech = pech_i, sech = sech_i,
        stringsAsFactors = FALSE
      )
    }
  }

  out = out[seq_len(k)]
  if (k == 0L) {
    return(data.frame(Period = character(), DMU = character(),
      mi = numeric(), ec = numeric(), tc = numeric(),
      pech = numeric(), sech = numeric(), stringsAsFactors = FALSE))
  }
  tibble::as_tibble(do.call(rbind, out))
}


#' Shephard distance function (raw distance, NOT Farrell efficiency).
#'
#' Output orientation: returns φ ≥ 1 (the Shephard output distance).
#' Input orientation:  returns θ ∈ (0,1] (the Shephard input distance).
#' @noRd
dist_shephard = function(x0, y0, X_ref, Y_ref, orientation, rts) {
  res = solve_radial_lp(X = X_ref, Y = Y_ref,
                        x0 = x0, y0 = y0,
                        orientation = orientation, rts = rts)
  if (orientation == "io") {
    if (is.na(res$theta)) return(NA_real_)
    res$theta
  } else {
    if (is.na(res$phi)) return(NA_real_)
    res$phi
  }
}
