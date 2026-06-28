# 评价模型 — 灰色预测
#
# 包含: GM11, GM1N, DGM21, verhulst, combine_preds
#

#' Grey Prediction Models
#'
#' @description
#' Implements grey prediction models for time series forecasting:
#' `GM11` applies the GM(1,1) model with level ratio test.
#' `GM1N` applies the GM(1,N) model with multiple related factors.
#' `DGM21` applies the DGM(2,1) model for second-order dynamics.
#' `verhulst` applies the Verhulst model for logistic growth.
#'
#' @param X For `GM11`, `DGM21`, `verhulst`: Numeric vector of original time series data.
#' @param dat For `GM1N`: Data frame or matrix, last column is characteristic series, others are related factors.
#' @param new_data For `GM1N`: Optional future values of related factors (1 row, m-1 columns).
#'
#' @return
#' For `GM11`: List with fitted values (`fitted`), next prediction (`pnext`), prediction function (`f`), matrix (`mat`), parameters (`u`), level ratios (`lambda`), and range (`rng`).
#' For `GM1N`: List with fitted values (`fitted`), posterior variance ratio (`C`), small error probability (`P`), and prediction function (`f`).
#' For `DGM21`, `verhulst`: List with fitted values (`fitted`), next prediction (`pnext`), prediction function (`f`), matrix (`mat`), and parameters (`u`).
#'
#' @examples
#' # Sample time series for GM11, DGM21, Verhulst
#' x = c(100, 120, 145, 175, 210)
#'
#' # GM11
#' result = GM11(x)
#' result$fitted    # Fitted values
#' result$pnext     # Next prediction
#' result$f(6:8)    # Predict next 3 periods
#'
#' # DGM21
#' x = c(2.874,3.278,3.39,3.679,3.77,3.8)
#' result = DGM21(x)
#' result$fitted    # Fitted values
#' result$pnext     # Next prediction
#' result$f(6:8)    # Predict next 3 periods
#'
#' # Verhulst
#' x = c(4.93,2.33,3.87,4.35,6.63,7.15,5.37,6.39,7.81,8.35)
#' result = verhulst(x)
#' result$fitted    # Fitted values
#' result$pnext     # Next prediction
#' result$f(6:8)    # Predict next 3 periods
#'
#' # Sample data for GM1N
#' data = data.frame(
#'   factor1 = c(50, 55, 60, 65, 70),
#'   factor2 = c(20, 22, 25, 28, 30),
#'   output = c(100, 120, 145, 175, 210)
#' )
#' result = GM1N(data)
#' result$fitted
#'
#' @name grey_models
NULL

# ============================================================================
# GM11 — GM(1,1) 灰色预测
# ============================================================================

#' @rdname grey_models
#' @export
GM11 = function(X) {
  if(!is.numeric(X) || is.matrix(X))
    stop("X must be a numeric vector.")
  if(length(X) < 4)
    stop("X must have length >= 4.")

  # Level ratio test
  n = length(X)
  lam = X[1:n-1] / X[2:n]
  rng = exp(c(-2,2) / (n + 1))
  if(all(rng[1]<lam & lam<rng[2])) {
    cat("Level ratio test passed!\n")
  } else {
    cat("Level ratio test failed!\n")
  }
  # GM(1,1) modeling
  X1 = cumsum(X)
  Z = (X1[1:n-1] + X1[2:n]) / 2
  B = cbind(-Z, rep(1,n-1))
  u = MASS::ginv(B) %*% X[2:n]
  k = 1:n
  X1h = (X[1]-u[2]/u[1])*exp(-u[1]*(k-1)) + u[2]/u[1]
  pred = c(X1h[1], diff(X1h))
  f = \(k) (X[1]-u[2]/u[1])*(1-exp(u[1]))*exp(-u[1]*(k-1))
  pnext = f(n+1)
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u,
       lambda = lam, rng = rng)
}


# ============================================================================
# GM1N — GM(1,N) 多因素灰色预测
# ============================================================================

#' @rdname grey_models
#' @export
GM1N = function(dat, new_data = NULL) {
  if(!is.data.frame(dat) && !is.matrix(dat))
    stop("dat must be a data frame or matrix.")

  X0 = as.matrix(dat)
  n = nrow(X0)
  m = ncol(X0)
  if(m < 2) stop("dat must have at least 2 columns.")
  if(n < 4) stop("dat must have at least 4 rows.")

  X = X0[, m]           # Characteristic series

  # Smoothness check
  X1 = cumsum(X)
  rho = X[2:n] / X1[1:(n-1)]
  rho_ratio = rho[2:(n-2)] / rho[1:(n-3)]
  smooth_flag = all(rho[3:(n-2)] <= 0.5 & rho_ratio < 1 & rho[n-1] <= 0.5)
  cat("Smoothness check:", ifelse(smooth_flag, "Passed", "Failed"), "\n")

  # Level ratio check
  lambds = X[1:(n-1)] / X[2:n]
  X_min = exp(-2 / (n + 1))
  X_max = exp(2 / (n + 1))
  level_flag = all(lambds >= X_min & lambds <= X_max)
  cat("Level ratio check:", ifelse(level_flag, "Passed", "Failed"), "\n")

  # Build GM(1,N) model
  X1 = apply(X0, 2, cumsum)
  Z = -0.5 * (X1[-n, m] + X1[-1, m])
  A = X0[-1, m]
  B = cbind(Z, X1[-1, -m])
  u = MASS::ginv(t(B) %*% B) %*% t(B) %*% A
  a = u[1]
  b = u[-1]

  # Prediction function
  f = function(k, X1_k) {
    (X[1] - (1 / a) * sum(X1_k * b)) * exp(-a * k) + (1 / a) * sum(X1_k * b)
  }

  # Fitted values
  X1h = sapply(1:n, function(k) f(k - 1, X1[k, -m]))
  fitted = c(X1h[1], diff(X1h))

  # Model evaluation
  S1 = sd(X)
  S2 = sd(X - fitted)
  C = S2 / S1
  Pe = mean(X - fitted)
  P = mean(abs(X - fitted - Pe) < 0.6745 * S1)

  # Predict next value
  pred = NULL
  if (!is.null(new_data)) {
    new_data = as.matrix(new_data)
    if (nrow(new_data) != 1 || ncol(new_data) != (m - 1)) {
      stop("new_data must be 1 row and ", m - 1, " columns matching the number of related factors.")
    }
    X1_last = X1[n, -m]
    X1_future = as.vector(new_data) + X1_last
    X1_ext = rbind(X1[, -m], X1_future)
    X1h_pred = sapply(0:n, function(k) f(k, X1_ext[k + 1, ]))
    pred = X1h_pred[n + 1] - X1h_pred[n]
  }

  list(fitted = fitted, pred = pred, vr = C, err = P, f = f)
}


# ============================================================================
# DGM21 — DGM(2,1) 二阶灰色预测
# ============================================================================

#' @rdname grey_models
#' @export
DGM21 = function(X) {
  if(!is.numeric(X) || is.matrix(X))
    stop("X must be a numeric vector.")
  if(length(X) < 4)
    stop("X must have length >= 4.")

  n = length(X)
  X1 = cumsum(X)
  AX1 = diff(X)
  B = cbind(-X[-1], rep(1, n-1))
  u = MASS::ginv(B) %*% AX1
  a = u[1]; b = u[2]
  x0 = X[1]
  f = \(k) (b/a^2 - x0/a) * exp(-a * (k-1)) + b/a * (k-1) +
    (1+a)/a * x0 - b/a^2
  X1h = f(1:(n+1))
  p = c(X1h[1], diff(X1h))
  pred = p[1:n]
  pnext = p[n+1]
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u)
}


# ============================================================================
# verhulst — Verhulst 灰色预测（Logistic 增长型）
# ============================================================================

#' @rdname grey_models
#' @export
#' @importFrom stats sd
verhulst = function(X) {
  if(!is.numeric(X) || is.matrix(X))
    stop("X must be a numeric vector.")
  if(length(X) < 4)
    stop("X must have length >= 4.")

  n = length(X)
  X1 = cumsum(X)
  Z = (X1[1:n-1] + X1[2:n]) / 2
  B = cbind(-Z, Z^2)
  u = MASS::ginv(B) %*% X[-1]
  f = \(k) u[1] * X[1] / (u[2] * X[1] + (u[1]-u[2]*X[1])*exp(u[1]*(k-1)))
  X1h = f(1:(n+1))
  p = c(X1h[1], diff(X1h))
  pred = p[1:n]
  pnext = p[n+1]
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u)
}


# ============================================================================
# combine_preds — 组合多模型预测
# ============================================================================

#' @title Combine Multiple Prediction Results
#'
#' @description Combines multiple prediction results (e.g., from grey prediction, time series, or machine learning models) into a single prediction using a similarity-based weighting approach, improving prediction accuracy.
#'
#' @param x Numeric vector, prediction results to be combined (length >= 2).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `a`: Numeric, the combined prediction value.
#'   \item `w`: Numeric vector, weights for each prediction in `x`, summing to 1.
#' }
#' @details
#' The function combines prediction results by constructing a similarity matrix based on cosine transformation of pairwise differences. Weights are derived from the principal eigenvector of the similarity matrix, ensuring predictions closer to each other have higher influence. For two predictions, equal weights (0.5, 0.5) are used. If all predictions are identical, equal weights are assigned.
#' Compatible with the `mathmodels` package for enhancing prediction models, including grey prediction, time series, or ensemble machine learning.
#' @examples
#' # Example: Combine three prediction results
#' preds = c(100, 102, 98)  # E.g., from grey prediction, ARIMA, or ML models
#' combine_preds(preds)

#' @export
#' @importFrom utils combn
combine_preds = function(x) {
  if(!is.numeric(x) || is.matrix(x))
    stop("x must be a numeric vector.")
  Lx = length(x)
  if (Lx < 2) stop("Input x must have length >= 2.")
  if(Lx == 2) {    # For two predictions, use equal weights
    a = mean(x)
    w = c(1/2,1/2)
  }
  IndCom = combn(1:Lx,2)
  d = abs(x[IndCom[1,]] - x[IndCom[2,]])
  maxd = max(d)
  idx = expand.grid(1:Lx,1:Lx)
  R = matrix(cos(pi*(x[idx[,1]]-x[idx[,2]]) / (2*maxd)),
             nrow = Lx, byrow = TRUE)
  rlt = eigen(R)
  w = rlt$vectors[,1] / sum(rlt$vectors[,1])
  a = sum(x * w)
  list(a = a, w = w)
}
