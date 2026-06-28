#' Markov Chain and Grey-Markov Prediction Models
#'
#' @description
#' Implements Markov chain prediction for state sequences and the Grey-Markov 
#' combined prediction model that integrates GM(1,1) with Markov chain 
#' correction.
#'
#' @param S For `markov_chain`: State sequence (vector or factor).
#' @param s0 For `markov_chain`: Initial state, must be one of the levels in `S`.
#' @param n_steps For `markov_chain`: Number of prediction steps (positive integer, default 5).
#' @param X For `GM11_markov`: Numeric vector of original time series data.
#' @param n_ahead For `GM11_markov`: Number of future periods to predict (default 3).
#' @param breaks For `GM11_markov`: Numeric vector of boundaries for state 
#'   classification based on relative error. Note: `-Inf` and `Inf` are 
#'   automatically prepended/appended internally.
#'
#' @return
#' \describe{
#'   \item{markov_chain}{Returns a list with:
#'     \itemize{
#'       \item \code{trans_mat}: Transition probability matrix.
#'       \item \code{pred_probs}: Matrix of state probabilities for each step.
#'       \item \code{pred_states}: Predicted states (most likely) for each step.
#'       \item \code{pi_final}: Ultimate stationary distribution.
#'     }
#'   }
#'   \item{GM11_markov}{Returns a data frame with columns:
#'     \itemize{
#'       \item \code{Period}: Time period labels (T1, T2, ..., T+n).
#'       \item \code{Raw}: Original values (NA for future periods).
#'       \item \code{GM11_fitted}: GM(1,1) fitted/predicted values.
#'       \item \code{err_state}: Relative error state (for historical) or 
#'              predicted state (for future).
#'       \item \code{adj_eff}: Adjustment effect (midpoint of state interval).
#'       \item \code{Markov_adj}: Markov chain adjusted values.
#'     }
#'   }
#' }
#'
#' @details
#' \code{markov_chain} constructs a transition probability matrix from a state 
#' sequence and performs multi-step prediction. For states that never appear 
#' as a starting state, equal transition probabilities are assigned to maintain 
#' row sums of 1. The ultimate stationary distribution is computed via 
#' eigenvalue decomposition.
#'
#' \code{GM11_markov} first fits a GM(1,1) model to the original series,
#' then classifies relative errors into states using the specified boundaries,
#' builds a Markov chain model on the error states, and corrects both historical
#' fitted values and future predictions.
#'
#' @examples
#' # --- Markov chain prediction ---
#' # Weather states: Rainy, Cloudy, Sunny
#' S = factor(c("Sunny", "Sunny", "Cloudy", "Rainy", "Sunny",
#'              "Cloudy", "Sunny", "Sunny", "Rainy", "Cloudy",
#'              "Sunny", "Cloudy", "Rainy", "Sunny", "Sunny",
#'              "Cloudy", "Sunny", "Rainy", "Cloudy", "Sunny"),
#'             levels = c("Rainy", "Cloudy", "Sunny"))
#' markov_chain(S, s0 = "Cloudy", n_steps = 3)

#'
#' # --- Grey-Markov prediction ---
#' X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
#' GM11_markov(X, n_ahead = 3)
#'
#' @name markov
NULL

#' @rdname markov
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


#' @rdname markov
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
