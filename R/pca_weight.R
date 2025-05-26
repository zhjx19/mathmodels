#' PCA-Based Weighting Method for Indicator Scoring
#'
#' Computes indicator weights using Principal Component Analysis (PCA).
#' The method extracts principal components and uses their variance contribution
#' to derive objective weights for indicators. Optionally handles positive/negative
#' directions of indicators.
#'
#' @param X A numeric data frame or matrix where rows represent samples and
#'          columns represent indicators.
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better).
#'              If not provided, all indicators are assumed to be positive.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of normalized weights for each indicator.}
#' \item{s}{Numeric vector of scores for each sample, scaled by 100.}
#' \item{lambda}{Eigenvalues (explained variance) of principal components.}
#' \item{B}{Loading matrix scaled by square root of eigenvalues.}
#' \item{beta}{Weight contributions from loadings and variance explained.}
#'
#' @importFrom psych principal
#' @export
#' @examples
#' # Example: Using PCA to compute indicator weights
#' ind = c("+","+","-","-")
#' pca_weight(iris[1:10, 1:4], ind, nfs = 2)


pca_weight = function(X, index = NULL, nfs = NULL) {
  # Compute indicator weights using Principal Component Analysis (PCA)
  # X: original data matrix, rows are samples, columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative
  # w: returned normalized weights, lambda: principal eigenvalues
  # B and beta: intermediate results (loadings matrix and weight contributions)

  if(is.null(index)) index = rep("+", ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], rescale, type = "-")

  # Perform PCA with all components and varimax rotation
  if(is.null(nfs)) nfs = ncol(X)
  pc = psych::principal(X, nfactors = nfs, rotate = "varimax")
  A = matrix(pc$loadings, ncol = pc$factors)
  lambda = pc$values[1:ncol(A)]

  # Scale loadings by sqrt of corresponding eigenvalues
  B = A / sqrt(matrix(rep(lambda, times = nrow(A)), ncol = ncol(A), byrow = TRUE))

  # Variance accounted by each component
  varP = pc$Vaccounted[2,]

  # Compute weighted contribution of each loading
  beta = abs(B %*% varP) / sum(varP)

  # Normalize weights
  w = beta / sum(beta)

  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)

  list(w = w[,1], s = s, lambda = lambda, B = B, beta = beta[,1])
}
#' @title PCA Weighting Method
