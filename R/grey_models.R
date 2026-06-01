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

#' @rdname grey_models
#' @export
GM11 = function(X) {
  # Implements GM(1,1) algorithm, input time series data
  # Returns fitted values (fitted), next prediction (pnext), prediction function (f)
  # Level ratios (lambda), acceptable range (rng)

  if(!is.numeric(X) || is.matrix(X))
    stop("X must be a numeric vector.")
  if(length(X) < 4)
    stop("X must have length >= 4.")

  # Level ratio test
  n = length(X)
  lam = X[1:n-1] / X[2:n]                 # Compute level ratios
  rng = exp(c(-2,2) / (n + 1))            # Acceptable range
  if(all(rng[1]<lam & lam<rng[2])) {
    cat("Level ratio test passed!\n")
  } else {
    cat("Level ratio test failed!\n")
  }
  # GM(1,1) modeling
  X1 = cumsum(X)                         # First-order accumulation
  Z = (X1[1:n-1] + X1[2:n]) / 2         # Mean of consecutive accumulated values
  B = cbind(-Z, rep(1,n-1))              # Construct matrix B
  u = MASS::ginv(B) %*% X[2:n]           # Compute parameters a, b
  k = 1:n
  # Predict accumulated values
  X1h = (X[1]-u[2]/u[1])*exp(-u[1]*(k-1)) + u[2]/u[1]
  pred = c(X1h[1], diff(X1h))            # Restore predicted values
  f = \(k) (X[1]-u[2]/u[1])*(1-exp(u[1]))*exp(-u[1]*(k-1))  # Prediction function
  pnext = f(n+1)
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u,
       lambda = lam, rng = rng)
}

#' @rdname grey_models
#' @export

GM1N = function(dat, new_data = NULL) {
  # Function: Simplified GM(1,N) grey model, predicts only one future value
  # Inputs:
  #   dat - data frame or matrix, last column is characteristic series, others are related factors
  #   new_data - future values of related factors (matrix or data frame, 1 row only)
  # Outputs:
  #   list containing fitted values (fitted), prediction (pred), variance ratio (vr), small error probability (err), and prediction function (f)

  if(!is.data.frame(dat) && !is.matrix(dat))
    stop("dat must be a data frame or matrix.")

  # Convert to matrix
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
  X1 = apply(X0, 2, cumsum)           # Accumulated Generating Operation (AGO)
  Z = -0.5 * (X1[-n, m] + X1[-1, m])  # Background values
  A = X0[-1, m]                       # Response variable
  B = cbind(Z, X1[-1, -m])            # Design matrix
  u = MASS::ginv(t(B) %*% B) %*% t(B) %*% A  # Parameter estimation
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

#' @rdname grey_models
#' @export
DGM21 = function(X) {
  # Implements DGM(2,1) algorithm, input time series data
  # Returns fitted values (fitted), next prediction (pnext), prediction function (f), matrix (mat), parameters (u)

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
  X1h = f(1:(n+1))           # Predict including next period
  p = c(X1h[1], diff(X1h))
  pred = p[1:n]               # Fitted values (first n)
  pnext = p[n+1]              # Next prediction
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u)
}

#' @rdname grey_models
#' @export
verhulst = function(X) {
  # Implements Verhulst model, input time series data
  # Returns fitted values (fitted), next prediction (pnext), prediction function (f), matrix (mat), parameters (u)

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
