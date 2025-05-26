#' Grey Prediction Models
#'
#' @description
#' Implements grey prediction models for time series forecasting:
#' \code{GM11} applies the GM(1,1) model with level ratio test.
#' \code{GM1N} applies the GM(1,N) model with multiple related factors.
#' \code{DGM21} applies the DGM(2,1) model for second-order dynamics.
#' \code{verhulst} applies the Verhulst model for logistic growth.
#'
#' @param X For \code{GM11}, \code{DGM21}, \code{verhulst}: Numeric vector of original time series data.
#' @param dat For \code{GM1N}: Data frame or matrix, last column is characteristic series, others are related factors.
#'
#' @return
#' For \code{GM11}: List with fitted values (\code{fitted}), next prediction (\code{pnext}), prediction function (\code{f}), matrix (\code{mat}), parameters (\code{u}), level ratios (\code{lambda}), and range (\code{rng}).
#' For \code{GM1N}: List with fitted values (\code{fitted}), posterior variance ratio (\code{C}), small error probability (\code{P}), and prediction function (\code{f}).
#' For \code{DGM21}, \code{verhulst}: List with fitted values (\code{fitted}), next prediction (\code{pnext}), prediction function (\code{f}), matrix (\code{mat}), and parameters (\code{u}).
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
GM1N = function(dat) {
  # Implements GM(1,N) algorithm
  # Input: data frame or matrix, last column is characteristic series, others are related factors
  # Returns fitted values (fitted), variance ratio (vr), error probability (err), prediction function (f)
  X0 = as.matrix(dat)
  n = nrow(X0)
  m = ncol(X0)
  X = X0[,m]
  # Level ratio test
  lam = X[1:n-1] / X[2:n]                # Compute level ratios
  rng = exp(c(-2,2) / (n + 1))           # Acceptable range
  if(all(rng[1]<lam & lam<rng[2])) {
    cat("Level ratio test passed!\n")
  } else {
    cat("Level ratio test failed!\n")
  }
  # GM(1,N) modeling
  X1 = apply(X0, 2, cumsum)
  Z = (X1[,m][-n] + X1[,m][-1]) / 2
  B = cbind(-Z, X1[-1,-m])
  Y = X[-1]
  u = MASS::ginv(B) %*% Y
  f = \(k, CX) {
    sbx = (1/u[1]) * sum(CX * u[2:m])
    (X[1] - sbx) * exp(-u[1] * (k-1)) + sbx
  }
  X1h = sapply(1:n, \(k) f(k, X1[k,-m]))
  pred = c(X1h[1], diff(X1h))       # Predicted values
  names(pred) = NULL
  # Model evaluation
  S1 = sd(X)
  S2 = sd(X - pred)
  C = S2 / S1           # Posterior variance ratio
  Pe = mean(X - pred)
  P = mean(abs((X - pred - Pe)) < 0.6745*S1)  # Small error probability
  list(fitted = pred, vr = C, err = P, f = f)
}

#' @rdname grey_models
#' @export
DGM21 = function(X) {
  # Implements DGM(2,1) algorithm, input time series data
  # Returns fitted values (pred), next prediction (pnext), prediction function (f), matrix (mat), parameters (u)
  n = length(X)
  X1 = cumsum(X)
  AX1 = diff(X)
  B = cbind(-X[-1], rep(1, n-1))
  u = MASS::ginv(B) %*% AX1
  a = u[1]; b = u[2]
  x0 = 2.874
  f = \(k) (b/a^2 - x0/a) * exp(-a * (k-1)) + b/a * (k-1) +
    (1+a)/a * x0 - b/a^2
  X1h = f(1:(n+1))           # Predict including next period
  p = c(X1h[1], diff(X1h))
  pred = p[-n]
  pnext = p[n]
  list(pred = fitted, pnext = pnext, f = f, mat = B, u = u)
}

#' @rdname grey_models
#' @export
verhulst = function(X) {
  # Implements Verhulst model, input time series data
  # Returns fitted values (pred), next prediction (pnext), prediction function (f), matrix (mat), parameters (u)
  n = length(X)
  X1 = cumsum(X)
  Z = (X1[1:n-1] + X1[2:n]) / 2
  B = cbind(-Z, Z^2)
  u = MASS::ginv(B) %*% X[-1]
  f = \(k) u[1] * X[1] / (u[2] * X[1] + (u[1]-u[2]*X[1])*exp(u[1]*(k-1)))
  X1h = f(1:(n+1))
  p = c(X1h[1], diff(X1h))
  pred = p[-n]
  pnext = p[n]
  list(fitted = pred, pnext = pnext, f = f, mat = B, u = u)
}
