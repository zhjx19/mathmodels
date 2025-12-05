#' TOPSIS Method for Multi-Criteria Decision Making
#'
#' @description
#' Implements the Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS)
#' to rank alternatives based on multiple criteria. The function normalizes the decision matrix using Min-Max method,
#' applies weights, and computes relative closeness to the ideal solution.
#'
#' @param X A numeric matrix or data frame where rows represent alternatives and columns represent criteria.
#' @param w A numeric vector of weights for each criterion. Must be non-negative and sum to 1.
#'          If not provided, equal weights are used.
#' @param index A character vector indicating the direction of each indicator:
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better).
#'              If `index = NULL` (default), all indicators are treated as `"+"`.
#'
#' @return A named numeric vector of relative closeness scores (in \eqn{[0, 1]}) for each alternative.
#'         Higher values indicate better alternatives.
#'         Names are taken from \code{rownames(X)} or default to "Sample1", "Sample2", etc.
#'
#' @details
#' The TOPSIS method ranks alternatives by:
#' \enumerate{
#'   \item Normalizing the decision matrix using Min-Max normalization.
#'   \item Applying weights to form a weighted normalized matrix.
#'   \item Identifying positive and negative ideal solutions based on indicator directions.
#'   \item Computing Euclidean distances to ideal solutions.
#'   \item Calculating relative closeness as \code{S0 / (S0 + Sstar)}, where \code{S0}
#'   is the distance to the negative ideal and \code{Sstar} is the distance to the positive ideal.
#' }
#' This implementation supports both positive and negative indicators via the \code{index} parameter.
#'
#' @export
#' @examples
#' A = data.frame(
#'   X1 = c(2, 5, 3),  # "+"
#'   X2 = c(8, 1, 6)   # "-"
#' )
#' w = c(0.6, 0.4)
#' idx = c("+","-")
#' topsis(A, w, idx)

topsis = function(X, w = NULL, index = NULL) {
  # Implements TOPSIS method
  # X: decision matrix
  # w: weight vector
  # index: direction of each indicator, "+" for positive, "-" for negative
  m = ncol(X)
  if(is.null(index)) index = rep("+", m)
  if(is.null(w)) w = rep(1/m, m)
  B = apply(X, 2, normalize)         # Normalize the decision matrix
  C = B %*% diag(w)                  # Weighted normalized matrix, or
  # C = t(t(B) * w)
  # C = B * matrix(rep(w, each = n), nrow = n)
  # sweep(B, 2, w, "*")
  # Calculate positive and negative ideal solutions
  max_vals = apply(C, 2, max)
  min_vals = apply(C, 2, min)
  Cstar = ifelse(index == "+", max_vals, min_vals)
  C0 = ifelse(index == "+", min_vals, max_vals)
  # Distance to positive ideal solution
  Sstar = apply(C, 1, \(x) norm(x - Cstar, "2"))
  # Distance to negative ideal solution
  S0 = apply(C, 1, \(x) norm(x - C0, "2"))
  S0 / (S0 + Sstar)                  # Relative closeness
}
