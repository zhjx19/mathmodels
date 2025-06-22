#' @title CRITIC Weight Method
#'
#' @description Computes objective weights of indicators and scores of samples using the CRITIC method.
#' The method considers both the contrast intensity (e.g., standard deviation or entropy)
#' and conflict among indicators (based on correlation) to determine indicator importance.
#' This version supports different methods for contrast intensity and correlation types.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#' @param index A character vector indicating the direction of each indicator:
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already normalized indicators (no rescaling will be applied).
#'
#'              If `index = NULL` (default), all indicators are treated as `NA`, meaning no normalization is performed.
#' @param method Character scalar; specifies the method used to compute contrast intensity.
#'               Options: \code{"std"} (standard deviation, default), or \code{"entropy"} (based on information redundancy).
#' @param cor_method Character scalar; specifies the method for computing correlations.
#'                   Options: \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#' @param epsilon A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors.
#'                Default is 0.002. Only used when \code{method = "entropy"}.

#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
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
#' critic_weight(X, index, method = "entropy")

critic_weight = function(X, index = NULL, method = "std",
                         cor_method = "pearson", epsilon = 0.002) {
  # Compute indicator weights and sample scores using the CRITIC method
  # X: original data matrix, rows are samples, columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative, NA means no rescaling.
  # method: "std" for standard deviation, "entropy" for entropy
  # cor_method: which correlation coefficient, "pearson"(default), "kendall" or "spearman"
  # epsilon: A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors, default is 0.002.
  # w: returned weights of indicators, s: returned scores of samples

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], rescale, type = "-")

  # Compute sum of conflicts (1 - correlation) for each indicator
  R = apply(1 - cor(X, method = cor_method), 2, sum)

  # Compute standard deviation or entropy of each indicator
  S = switch (method,
    "std" = apply(X, 2, sd),
    "entropy" = {
      X[X == 0] = epsilon
      X[X == 1] = 1 - epsilon
      P = data.frame(lapply(X, \(x) x / sum(x)))
      e = sapply(P, function(x) - sum(x * log(x)) /log(nrow(P)))
      1 - e
    })

  # Compute indicator importance (C)
  C = S * R
  # Compute weights
  w = C / sum(C)
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)

  list(w = w, s = s)
}
#' @title CRITIC Method
