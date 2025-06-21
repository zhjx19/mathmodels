#' @title CRITIC Weight Method
#'
#' @description Computes objective weights of indicators and scores of samples using the CRITIC method.
#' The method considers both the variance (contrast intensity) and correlation
#' (conflict among indicators) to determine indicator importance.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#'
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already normalized indicators (no rescaling will be applied).
#'
#'              If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no normalization or rescaling is performed.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}#'
#'
#' @export
#' @examples
#' # Example: Using CRITIC method on a simple dataset
#' X = data.frame(
#'   x1 = c(3, 5, 2, 7),
#'   x2 = c(10, 20, 15, 25)
#' )
#' index = c("+", "-")
#' critic_weight(X, index)

critic_weight = function(X, index = NULL) {
  # Compute indicator weights and sample scores using the CRITIC method
  # X: original data matrix, rows are samples, columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative, NA means no rescaling.
  # w: returned weights of indicators, s: returned scores of samples

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], rescale, type = "-")

  # Compute standard deviation of each indicator
  S = apply(X, 2, sd)

  # Compute sum of conflicts (1 - correlation) for each indicator
  R = apply(1 - cor(X), 2, sum)

  # Compute indicator importance (C)
  C = S * R

  # Compute weights
  w = C / sum(C)

  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)

  list(w = w, s = s)
}
#' @title CRITIC Method
