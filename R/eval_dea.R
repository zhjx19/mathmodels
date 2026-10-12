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
#' Returns matrices in m\eqn{\times}n format (inputs/outputs as rows, DMUs as columns),
#' matching the MATLAB convention used in DEA linear programming.
#' @noRd
extract_dea_data = function(data, inputs, outputs, ud_outputs = NULL) {
  .check_matrix_like(data, "`data` must be a data frame.", df_only = TRUE)
  if (ncol(data) < 3) stop("`data` must have at least 3 columns (DMU + >=1 input + >=1 output).")

  dmu_names = data[[1]]
  X = t(as.matrix(data[, inputs, drop = FALSE]))      # m × n
  Y_all = t(as.matrix(data[, outputs, drop = FALSE]))   # s × n

  if (anyNA(X) || anyNA(Y_all)) stop("Input/output columns must not contain NA values.")
  if (any(X <= 0) || any(Y_all <= 0))
    stop("DEA requires strictly positive inputs and outputs.", call. = FALSE)

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
#' @param X Input matrix (m \eqn{\times} n_ref)
#' @param Y Desirable output matrix (s1 \eqn{\times} n_ref)
#' @param x0 Target DMU input vector (length m)
#' @param y0 Target DMU desirable output vector (length s1)
#' @param orientation "io" or "oo"
#' @param rts "crs" or "vrs"
#' @param Y_b Undesirable output matrix (s2 \eqn{\times} n_ref), NULL if none.
#'   Treated as an input (Hailu & Veeman 2001): io — bad outputs contract with
#'   \eqn{\theta}; oo — bad outputs held fixed while good outputs expand.
#' @param y0_b Target DMU undesirable output vector (length s2), NULL if none.
#' @return list(theta, phi, lambda, s_in, s_out, s_b)
#'   theta = Farrell efficiency (0~1); phi = Shephard output distance (\eqn{\ge 1}, NA for io)
#' @noRd
solve_radial_lp = function(X, Y, x0, y0, orientation = "io", rts = "vrs",
                           Y_b = NULL, y0_b = NULL) {
  n_ref = ncol(X)
  m = length(x0)
  q = length(y0)
  s2 = if (is.null(Y_b)) 0L else nrow(Y_b)

  lp = lpSolveAPI::make.lp(0, n_ref + 1)

  if (orientation == "io") {
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), 1))  # min θ
    for (k in seq_len(m))
      lpSolveAPI::add.constraint(lp, c(X[k, ], -x0[k]), "<=", 0)
    for (k in seq_len(q))
      lpSolveAPI::add.constraint(lp, c(-Y[k, ], 0), "<=", -y0[k])
    for (k in seq_len(s2))
      lpSolveAPI::add.constraint(lp, c(Y_b[k, ], -y0_b[k]), "<=", 0)
  } else {
    lpSolveAPI::set.objfn(lp, c(rep(0, n_ref), -1))  # min -φ ≡ max φ
    for (k in seq_len(m))
      lpSolveAPI::add.constraint(lp, c(X[k, ], 0), "<=", x0[k])
    for (k in seq_len(q))
      lpSolveAPI::add.constraint(lp, c(-Y[k, ], y0[k]), "<=", 0)
    for (k in seq_len(s2))
      lpSolveAPI::add.constraint(lp, c(Y_b[k, ], 0), "<=", y0_b[k])
  }

  if (rts == "vrs")
    lpSolveAPI::add.constraint(lp, c(rep(1, n_ref), 0), "=", 1)

  lpSolveAPI::set.bounds(lp, lower = rep(0, n_ref + 1))

  status = lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(theta = NA_real_, phi = NA_real_,
                lambda = rep(NA_real_, n_ref),
                s_in = rep(NA_real_, m), s_out = rep(NA_real_, q),
                s_b = rep(NA_real_, s2)))
  }

  sol = lpSolveAPI::get.variables(lp)
  lambda = sol[1:n_ref]
  theta_raw = sol[n_ref + 1]

  if (orientation == "io") {
    s_in  = pmax(0, theta_raw * x0 - as.vector(X %*% lambda))
    s_out = pmax(0, as.vector(Y %*% lambda) - y0)
    s_b   = if (s2 > 0) pmax(0, theta_raw * y0_b - as.vector(Y_b %*% lambda)) else numeric(0)
    eff = theta_raw
    phi = NA_real_
  } else {
    phi = theta_raw
    s_in  = pmax(0, x0 - as.vector(X %*% lambda))
    s_out = pmax(0, as.vector(Y %*% lambda) - phi * y0)
    s_b   = if (s2 > 0) pmax(0, y0_b - as.vector(Y_b %*% lambda)) else numeric(0)
    eff = 1 / phi
  }

  list(theta = eff, phi = phi,
       lambda = lambda, s_in = s_in, s_out = s_out, s_b = s_b)
}

#' Solve SBM LP (Tone 2001, Charnes-Cooper transformed) for a single DMU
#'
#' @param X Input matrix (m \eqn{\times} n_ref)
#' @param Y_g Desirable output matrix (s1 \eqn{\times} n_ref)
#' @param Y_b Undesirable output matrix (s2 \eqn{\times} n_ref), NULL if none
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
  # Tone (2004) normalization weights: desirable slacks 1/s1, undesirable 1/s2.
  # With no undesirable outputs (s2 = 0) this reduces to the plain SBM 1/s1.
  c_g = 1 / s1
  c_b = if (s2 > 0) 1 / s2 else 0

  lp = lpSolveAPI::make.lp(0, nvar)

  obj = numeric(nvar)
  obj[1:m] = -1 / (m * x0)
  obj[nvar] = 1
  lpSolveAPI::set.objfn(lp, obj)

  # Constraint 1: output slack weighted sum + t = 1
  con1 = numeric(nvar)
  con1[(m + 1):(m + s1)] = c_g / y0_g
  if (s2 > 0) con1[(m + s1 + 1):(m + s1 + s2)] = c_b / y0_b
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

  targets_out = dd$Y_g + t(s_out_mat[, seq_len(dd$s1), drop = FALSE])
  if (dd$s2 > 0) {
    targets_out = rbind(
      targets_out,
      dd$Y_b - t(s_out_mat[, (dd$s1 + 1):(dd$s1 + dd$s2), drop = FALSE])
    )
  }

  list(
    efficiencies = theta_vec,
    lambdas = lambda_mat,
    slacks = list(
      inputs   = s_in_mat,
      outputs  = s_out_mat
    ),
    targets = list(
      inputs   = dd$X - t(s_in_mat),
      outputs  = targets_out
    ),
    returns = structure(rep(rts_label, dd$n), names = dd$dmu),
    model = "DEA",
    orientation = rlang::`%||%`(attr(dd, "orientation"), "io"),
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
#'   For output orientation: value = 1/\eqn{\phi}, where \eqn{\phi \ge 1} is the Shephard distance.
#' @param rts `"vrs"` (variable returns to scale, default) or `"crs"` (constant).
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

  failed = character(0)
  for (i in seq_len(n)) {
    y_g = dd$Y_g[, i]
    y_b = if (s2 > 0) dd$Y_b[, i] else NULL
    res = solve_radial_lp(X = dd$X, Y = dd$Y_g,
                          x0 = dd$X[, i], y0 = y_g,
                          orientation = orientation, rts = rts,
                          Y_b = dd$Y_b, y0_b = y_b)
    if (!is.na(res$theta)) {
      theta_vec[i] = res$theta
      lambda_mat[i, ] = res$lambda
      s_in_mat[i, ] = res$s_in
      s_out_mat[i, seq_len(s1)] = res$s_out
      if (s2 > 0) s_out_mat[i, (s1 + 1):(s1 + s2)] = res$s_b
    } else {
      failed = c(failed, dd$dmu[i])
    }
  }
  if (length(failed))
    warning("basic_DEA: LP infeasible for DMU(s): ",
            paste(failed, collapse = ", "), "; efficiency set to NA.",
            call. = FALSE)

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
  if (n <= 1)
    warning("super_DEA needs at least 2 DMUs; the result is NA.", call. = FALSE)

  theta_vec = setNames(rep(NA_real_, n), dd$dmu)
  lambda_mat = matrix(NA_real_, n, n,
    dimnames = list(dd$dmu, dd$dmu))
  s_in_mat  = matrix(NA_real_, n, m,
    dimnames = list(dd$dmu, colnames(t(dd$X))))
  s_out_mat = matrix(NA_real_, n, s1,
    dimnames = list(dd$dmu, colnames(t(dd$Y_g))))

  failed = character(0)
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
    } else {
      failed = c(failed, dd$dmu[i])
    }
  }
  if (length(failed))
    warning("super_DEA: LP infeasible for DMU(s): ",
            paste(failed, collapse = ", "), "; efficiency set to NA.",
            call. = FALSE)

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
  if (!identical(orientation, "io"))
    stop("SBM models are non-oriented; `orientation` is not applicable.",
         call. = FALSE)
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

#' Solve super-efficiency SBM LP (Tone 2002, Charnes-Cooper transformed)
#' for a single DMU with the DMU excluded from the reference set.
#'
#' LP: min tau = t + (1/m) sum_i S_i^- / x_i0
#'     s.t. t x_i0 + S_i^- >= sum_{j!=0} lambda_j x_ij      (inputs)
#'          t y_r0 - S_r^+ <= sum_{j!=0} lambda_j y_rj      (desirable outputs)
#'          t - (1/s1) sum_r S_r^+ / y_r0 = 1               (normalization)
#'          VRS: sum_{j!=0} lambda_j = 1
#'
#' Returns tau >= 1 for SBM-efficient DMUs and exactly 1 for SBM-inefficient
#' ones (a dominating combination always exists in the remaining set).
#'
#' @param X Input matrix (m \eqn{\times} n_ref), DMU0 column already removed
#' @param Y_g Desirable output matrix (s1 \eqn{\times} n_ref), DMU0 removed
#' @param x0 Target DMU input vector (length m)
#' @param y0_g Target DMU desirable output vector (length s1)
#' @param rts "crs" or "vrs"
#' @return list(rho, s_in, s_g, lambda) with slacks/lambda in original units
#' @noRd
solve_super_sbm_lp = function(X, Y_g, x0, y0_g, rts = "vrs") {
  n_ref = ncol(X)
  m = length(x0)
  s1 = length(y0_g)

  nvar = n_ref + m + s1 + 1
  ii_lambda = 1:n_ref
  ii_sin = n_ref + (1:m)
  ii_sout = n_ref + m + (1:s1)
  ii_t = nvar

  lp = lpSolveAPI::make.lp(0, nvar)

  obj = numeric(nvar)
  obj[ii_sin] = 1 / (m * x0)
  obj[ii_t] = 1
  lpSolveAPI::set.objfn(lp, obj)

  # Normalization: t - (1/s1) * sum_r S_r^+ / y_r0 = 1
  row = numeric(nvar)
  row[ii_sout] = -1 / (s1 * y0_g)
  row[ii_t] = 1
  lpSolveAPI::add.constraint(lp, row, "=", 1)

  # Input rows: t * x_i0 + S_i^- >= sum_j lambda_j x_ij
  for (k in seq_len(m)) {
    row = numeric(nvar)
    row[ii_lambda] = -X[k, ]
    row[ii_sin[k]] = 1
    row[ii_t] = x0[k]
    lpSolveAPI::add.constraint(lp, row, ">=", 0)
  }

  # Desirable output rows: t * y_r0 - S_r^+ <= sum_j lambda_j y_rj
  for (k in seq_len(s1)) {
    row = numeric(nvar)
    row[ii_lambda] = Y_g[k, ]
    row[ii_sout[k]] = -1
    row[ii_t] = -y0_g[k]
    lpSolveAPI::add.constraint(lp, row, ">=", 0)
  }

  # VRS: sum lambda = 1
  if (rts == "vrs") {
    row = numeric(nvar)
    row[ii_lambda] = 1
    lpSolveAPI::add.constraint(lp, row, "=", 1)
  }

  lpSolveAPI::set.bounds(lp, lower = rep(0, nvar))

  status = lpSolveAPI::solve.lpExtPtr(lp)
  if (status != 0) {
    return(list(rho = NA_real_, s_in = rep(NA_real_, m),
                s_g = rep(NA_real_, s1), lambda = rep(NA_real_, n_ref)))
  }

  sol = lpSolveAPI::get.variables(lp)
  t_val = sol[ii_t]
  rho = lpSolveAPI::get.objective(lp)

  list(rho = rho,
       s_in = sol[ii_sin] / t_val,
       s_g = sol[ii_sout] / t_val,
       lambda = sol[ii_lambda] / t_val)
}

#' @describeIn DEA Super-efficiency SBM (Tone 2002; self excluded, no
#'   undesirable outputs)
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
  if (!identical(orientation, "io"))
    stop("SBM models are non-oriented; `orientation` is not applicable.",
         call. = FALSE)
  rts = match.arg(rts, c("vrs", "crs"))

  dd = extract_dea_data(data, inputs, outputs, ud_outputs = NULL)
  n = dd$n; m = dd$m; s1 = dd$s1
  if (n <= 1)
    warning("super_SBM needs at least 2 DMUs; the result is NA.", call. = FALSE)

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
    res = solve_super_sbm_lp(
      X = dd$X[, idx_ref, drop = FALSE],
      Y_g = dd$Y_g[, idx_ref, drop = FALSE],
      x0 = dd$X[, i], y0_g = dd$Y_g[, i],
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
#' @param period Column name or index identifying the time period variable.
#'   Numeric periods are sorted numerically; character periods that all parse
#'   as numbers are sorted numerically as well.
#' @param type1 `"cont"` (contemporaneous), `"seq"` (sequential), or `"glob"`
#'   (global, default). With `"seq"`, the old-technology term of the
#'   technical-change component at transition t -> t+1 is evaluated against
#'   the pooled frontier of periods 1..t (sequential frontier); no scale
#'   decomposition is defined (`sech = NA`).
#' @param type2 `"fgnz"` (Färe-Grosskopf-Norris-Zhang) or `"rd"` (default).
#'   `"rd"` reports a scale-efficiency component in the style of the **deaR**
#'   implementation (closest to Ray-Desli 1997: the scale term is measured at
#'   the period-t bundle relative to the mixed-period bundle); it is not the
#'   symmetric Ray-Desli decomposition, and `mi` is identical to the `"fgnz"`
#'   value.
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

  .check_matrix_like(data, "`data` must be a data frame.", df_only = TRUE)
  if (ncol(data) < 4)       stop("data must have >= 4 columns.")

  dmu_vec    = data[[1]]
  period_raw = data[[period]]
  X_all = t(as.matrix(data[, inputs,  drop = FALSE]))
  Y_all = t(as.matrix(data[, outputs, drop = FALSE]))
  if (anyNA(X_all) || anyNA(Y_all))
    stop("Input/output columns must not contain NA values.", call. = FALSE)
  if (any(X_all <= 0) || any(Y_all <= 0))
    stop("DEA requires strictly positive inputs and outputs.", call. = FALSE)

  # Sort periods numerically when the column is character but numeric-looking
  # (lexicographic sorting would put "10" before "2" and invert the index).
  if (is.factor(period_raw)) period_raw = as.character(period_raw)
  if (is.character(period_raw)) {
    pv = suppressWarnings(as.numeric(period_raw))
    if (anyNA(pv)) {
      warning("period values are non-numeric; sorting lexicographically.",
              call. = FALSE)
    } else {
      period_raw = pv
    }
  }
  period_vec = period_raw

  periods = sort(unique(period_vec))
  nt = length(periods)
  if (nt < 2) stop("At least 2 periods required.")

  dmu0 = dmu_vec[period_vec == periods[1]]
  nd   = length(dmu0)

  # Panel hygiene warnings
  for (p in periods) {
    if (anyDuplicated(dmu_vec[period_vec == p]))
      warning("Duplicate DMU names in period ", p, "; first occurrence used.",
              call. = FALSE)
  }
  missing_dmu = sum(vapply(periods[-1], function(p) {
    sum(!dmu0 %in% dmu_vec[period_vec == p])
  }, numeric(1)))
  if (missing_dmu > 0)
    warning(missing_dmu, " DMU-period combination(s) missing; ",
            "affected indices are NA.", call. = FALSE)

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

  # ---- 3b. Sequential distances (type1 = "seq") ----
  # Reference technology at transition t -> t+1 pools all DMUs of periods
  # 1..t (Grosskopf-style sequential frontier). Only the "old technology"
  # side of the technical-change term becomes cumulative; the contemporaneous
  # terms stay period-specific, so no self-reference occurs.
  D12cum_crs = matrix(NA_real_, nt - 1, nd, dimnames = list(periods[-1], dmu0))
  if (type1 == "seq") {
    for (t_idx in seq_len(nt - 1)) {
      idx_cum = which(period_vec %in% periods[seq_len(t_idx)])
      X_cum = X_all[, idx_cum, drop = FALSE]
      Y_cum = Y_all[, idx_cum, drop = FALSE]
      idx_t1 = which(period_vec == periods[t_idx + 1])
      dmu_t1 = dmu_vec[idx_t1]
      for (j in seq_len(nd)) {
        ii_t1 = match(dmu0[j], dmu_t1)
        if (is.na(ii_t1)) next
        D12cum_crs[t_idx, j] = dist_shephard(
          X_all[, idx_t1[ii_t1]], Y_all[, idx_t1[ii_t1]],
          X_cum, Y_cum, orientation, "crs")
      }
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
  E12cum = to_farrell(D12cum_crs)

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

      } else if (type1 == "seq") {
        # Sequential: TC uses the cumulative frontier of periods 1..t for the
        # old-technology evaluation; EC stays contemporaneous. No scale
        # decomposition is defined for the sequential index.
        tc_i = (E12cum[t_idx, j] / E[t_idx + 1, j]) *
               (E[t_idx, j] / E21[t_idx, j])
        ec_i = E[t_idx + 1, j] / E[t_idx, j]
        mi_i = ec_i * tc_i
        if (rts == "crs") {
          pech_i = NA_real_; sech_i = NA_real_
        } else {
          pech_i = Ev[t_idx + 1, j] / Ev[t_idx, j]  # reference only
          sech_i = NA_real_
        }

      } else {
        # cont: Färe et al. (1994) geometric mean
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
#' Output orientation: returns \eqn{\phi \ge 1} (the Shephard output distance).
#' Input orientation:  returns \eqn{\theta \in (0,1]} (the Shephard input distance).
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
