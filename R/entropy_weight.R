#' Entropy Weight Method for Indicator Weighting and Sample Scoring
#'
#' Computes the weights of indicators and scores of samples based on the entropy method.
#' This method objectively determines the importance of each indicator according to
#' the amount of information it contains.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better) and `"-"` for negative
#'              indicators (lower is better). If not provided, all indicators are assumed
#'              to be positive.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
#'
#' @export
#' @examples
#' # Example: Using entropy weight method on a simple dataset
#' X = data.frame(
#'   x1 = c(3, 5, 2, 7),
#'   x2 = c(10, 20, 15, 25)
#' )
#' index = c("+", "-")
#' entropy_weight(X, index)

entropy_weight = function(X, index = NULL) {
  # Compute weights of indicators and scores of samples using entropy method
  # X: original data matrix, rows are samples and columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative
  # w: returned weights of indicators, s: returned scores of samples
  if(is.null(index)) index = rep("+", ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale, a = 0.002, b = 0.996)
  X[,neg] = lapply(X[,neg, drop = FALSE], rescale, type = "-", a = 0.002, b = 0.996)
  # Compute proportion p(i,j) of sample i in indicator j
  P = data.frame(lapply(X, \(x) x / sum(x)))
  # Compute entropy e(j) for each indicator j
  e = sapply(P, \(x) sum(x * log(x)) *(-1/log(nrow(P))))
  d = 1 - e          # Compute redundancy degree
  w = d / sum(d)     # Compute weight vector
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)
  list(w = w, s = s)
}
#' @title Entropy Weight Method
