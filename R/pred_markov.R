#' Markov Chain Prediction
#'
#' @description
#' Constructs a transition probability matrix from a state sequence and performs
#' multi-step Markov chain prediction.
#'
#' @param S State sequence, supplied as a vector or factor.
#' @param s0 Initial state, must be one of the state levels in `S`.
#' @param n_steps Number of prediction steps. Must be a positive integer.
#'
#' @return A list with components:
#' \describe{
#'   \item{trans_mat}{Transition probability matrix.}
#'   \item{pred_probs}{Matrix of state probabilities for each prediction step.}
#'   \item{pred_states}{Predicted states, taking the most likely state at each step.}
#'   \item{pi_final}{Ultimate stationary distribution.}
#' }
#'
#' @details
#' For states that never appear as a starting state, equal transition
#' probabilities are assigned to keep each row sum equal to 1. The ultimate
#' stationary distribution is computed using eigenvalue decomposition.
#'
#' @examples
#' # Weather states: Rainy, Cloudy, Sunny
#' S = factor(c("Sunny", "Sunny", "Cloudy", "Rainy", "Sunny",
#'              "Cloudy", "Sunny", "Sunny", "Rainy", "Cloudy",
#'              "Sunny", "Cloudy", "Rainy", "Sunny", "Sunny",
#'              "Cloudy", "Sunny", "Rainy", "Cloudy", "Sunny"),
#'             levels = c("Rainy", "Cloudy", "Sunny"))
#' markov_chain(S, s0 = "Cloudy", n_steps = 3)
#'
#' @export
markov_chain = function(S, s0, n_steps = 5) {
  # Input validation
  if(!is.vector(S) && !is.factor(S))
    stop("S must be a vector or factor.")
  if(length(S) < 2)
    stop("S must have length >= 2.")
  
  if(!is.factor(S)) {
    S = factor(S, levels = sort(unique(S)))
  }
  lvl = levels(S)
  
  if(!is.character(s0) || length(s0) != 1)
    stop("s0 must be a single character string.")
  if(!s0 %in% lvl)
    stop("s0 must be one of the state levels: ", 
         paste(lvl, collapse = ", "), ".")
  
  if(!is.numeric(n_steps) || length(n_steps) != 1 || 
     n_steps < 1 || n_steps != round(n_steps))
    stop("n_steps must be a positive integer.")
  
  # Build transition matrix (table() with factor ensures square matrix 
  # covering all levels)
  freq = table(S[-length(S)], S[-1])
  P = prop.table(freq, 1)
  # For states that only appear as destination but never as source:
  # assign equal probabilities to keep row sums = 1
  P[is.nan(P)] = 1 / ncol(P)
  
  # Iterative prediction
  pred = matrix(NA_real_, nrow = n_steps + 1, ncol = length(lvl))
  pred[1, ] = as.numeric(lvl == s0)
  
  for(i in 1:n_steps) {
    pred[i + 1, ] = pred[i, ] %*% P
  }
  pred = pred[-1, , drop = FALSE]
  rownames(pred) = paste0("T", 1:n_steps)
  colnames(pred) = lvl
  
  # Most likely state at each step
  pred_states = lvl[max.col(pred, ties.method = "first")]
  
  # Ultimate stationary distribution via eigenvalue decomposition
  eig = eigen(t(P))
  idx = which.max(Re(eig$values))
  v = eig$vectors[, idx]
  pi_final = abs(Re(v)) / sum(abs(Re(v)))
  names(pi_final) = lvl
  
  list(trans_mat = P, pred_probs = pred,
       pred_states = pred_states, pi_final = pi_final)
}


#' Grey-Markov Prediction Model
#'
#' @description
#' Fits a GM(1,1) model and then uses a Markov chain on relative error states
#' to correct both historical fitted values and future predictions.
#'
#' @param X Numeric vector of original time series data.
#' @param n_ahead Number of future periods to predict. Must be a positive
#'   integer.
#' @param breaks Numeric vector of boundaries for classifying relative errors
#'   into states. `-Inf` and `Inf` are prepended/appended internally.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Period}{Time period labels, such as `T1`, `T2`, ..., `T+n`.}
#'   \item{Raw}{Original values; future periods are `NA`.}
#'   \item{GM11_fitted}{GM(1,1) fitted or predicted values.}
#'   \item{err_state}{Relative error state for historical periods, and predicted
#'     state for future periods.}
#'   \item{adj_eff}{Adjustment effect, defined as the midpoint of the state
#'     interval.}
#'   \item{Markov_adj}{Markov-chain-adjusted fitted or predicted values.}
#' }
#'
#' @details
#' The function first fits a GM(1,1) model to the original series, classifies
#' relative errors into states using `breaks`, builds a Markov chain on these
#' error states, and finally adjusts the GM(1,1) fitted and predicted values.
#'
#' @examples
#' X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
#' GM11_markov(X, n_ahead = 3)
#'
#' @export
GM11_markov = function(X, n_ahead = 3,
                       breaks = c(-0.10, -0.05, -0.02, 0, 0.02, 0.05, 0.10)) {
  # Input validation
  if(!is.numeric(X) || is.matrix(X))
    stop("X must be a numeric vector.")
  n = length(X)
  if(n < 4)
    stop("GM(1,1) requires at least 4 data points.")
  
  if(!is.numeric(n_ahead) || length(n_ahead) != 1 ||
     n_ahead < 1 || n_ahead != round(n_ahead))
    stop("n_ahead must be a positive integer.")
  
  if(!is.numeric(breaks) || length(breaks) < 1)
    stop("breaks must be a numeric vector.")
  
  # 1. Grey prediction GM(1,1)
  rlt = GM11(X)
  pred_hist = rlt$fitted
  pred_fut = rlt$f(n + 1:n_ahead)
  
  # 2. Relative error and state classification
  err = (X - pred_hist) / X
  
  mids = c(breaks[1], (breaks[-length(breaks)] + breaks[-1]) / 2, 
           breaks[length(breaks)])
  
  K = length(breaks) + 1
  S = cut(err, breaks = c(-Inf, breaks, Inf), labels = 1:K, 
          right = FALSE, include.lowest = TRUE) |>
    factor(levels = 1:K)
  
  # 3. Markov correction for historical fitted values (single-step rolling)
  # Period 1 cannot be corrected; periods 2~n use t-1's true state to 
  # predict t's state for correction
  ps_fit = sapply(as.character(S[-n]), 
                  \(x) markov_chain(S, x, 1)$pred_states) |> 
    as.numeric()
  
  fitted_new = c(pred_hist[1], pred_hist[-1] / (1 - mids[ps_fit]))
  
  # 4. Markov correction for future predictions
  last_state = as.character(S[n])
  mc_fut = markov_chain(S, last_state, n_ahead)
  ps_fut = as.numeric(mc_fut$pred_states)
  pnext_new = pred_fut / (1 - mids[ps_fut])
  
  # 5. Assemble result data frame
  data.frame(
    Period = c(paste0("T", 1:n), paste0("T+", 1:n_ahead)),
    Raw = c(X, rep(NA, n_ahead)),
    GM11_fitted = round(c(pred_hist, pred_fut), 2),
    err_state = c(as.character(S), mc_fut$pred_states),
    adj_eff = round(c(mids[as.numeric(S)], mids[ps_fut]), 4),
    Markov_adj = round(c(fitted_new, pnext_new), 2))
}
