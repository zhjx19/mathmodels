#' @title Entropy Weight Method
#'
#' @description Computes the weights of indicators and scores of samples based on the entropy method.
#' This method objectively determines the importance of each indicator according to
#' the amount of information it contains.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#'
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already normalized indicators (no rescaling will be applied,
#'              but minor adjustments will still be made to avoid log(0) errors).
#'              If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no normalization or rescaling is performed,
#'              but a small adjustment is still applied to prevent log(0) errors.
#' @param epsilon A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors.
#'                Default is 0.002.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
#'
#' @export
#'
#' @examples
#' X = data.frame(
#'   x1 = c(3, 5, 2, 7),
#'   x2 = c(10, 20, 15, 25)
#' )
#' index = c("+", "-")
#' entropy_weight(X, index)
#'

entropy_weight = function(X, index = NULL, epsilon = 0.002) {
  # Compute weights of indicators and scores of samples using entropy method
  # X: original data matrix, rows are samples and columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative, NA means no rescaling.
  # w: returned weights of indicators, s: returned scores of samples
  if(is.null(index)) index = rep(NA, ncol(X))
  if(any(is.na(index))) {
    X_na = X[is.na(index)]
    X_na[X_na == 0] = epsilon
    X_na[X_na == 1] = 1- epsilon
    X[is.na(index)] = X_na
  }
  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE],
                   function(x) rescale(x, a = epsilon, b = 1-epsilon))
  X[,neg] = lapply(X[,neg, drop = FALSE],
                   function(x) rescale(x, type = "-", a = epsilon, b = 1-epsilon))
  # Compute proportion p(i,j) of sample i in indicator j
  P = data.frame(lapply(X, function(x) x / sum(x)))
  # Compute entropy e(j) for each indicator j
  e = sapply(P, \(x) -sum(x * log(x)) / log(nrow(P)))
  d = 1 - e          # Compute redundancy degree
  w = d / sum(d)     # Compute weight vector
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)
  list(w = w, s = s)
}
#' @title Entropy Weight Method
