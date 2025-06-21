#' @title PCA-Based Weighting Method
#'
#' @description Computes indicator weights using Principal Component Analysis (PCA).
#' The method extracts principal components and uses their variance contribution
#' to derive objective weights for indicators. Optionally handles positive/negative
#' directions of indicators, and supports pre-standardized data.
#'
#' @param X A numeric data frame or matrix where rows represent samples and
#'          columns represent indicators.
#'
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already standardized indicators (no standardization will be applied).
#'
#'              If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no standardization is performed.
#'
#' @param nfs Number of principal components to use; by default, all are used.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of normalized weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample, scaled by 100.}
#'
#' \item{lambda}{Eigenvalues of principal components (explained variance).}
#'
#' \item{B}{Loading matrix scaled by square root of eigenvalues.}
#'
#' \item{beta}{Weight contributions derived from loadings and variance explained.}
#'
#' @export
#' @examples
#' # Example: Using PCA to compute indicator weights
#' ind = c("+","+","-","-")
#' pca_weight(iris[1:10, 1:4], ind, nfs = 2)

pca_weight = function(X, index = NULL, nfs = NULL) {
  # Compute indicator weights using Principal Component Analysis (PCA)
  # X: original data matrix, rows are samples, columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative
  # w: returned normalized weights, lambda: eigenvalues
  # B and beta: intermediate results (loadings matrix and weight contributions)

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Standardize data
  X = scale(X)
  X[,neg] = -X[,neg]

  # Perform PCA with all components and varimax rotation
  if(is.null(nfs)) nfs = ncol(X)

  # Replace psych::principal with prcomp + varimax
  pc = prcomp(X, center = FALSE, scale. = FALSE) # Already scaled above
  loadings = pc$rotation[, 1:nfs, drop = FALSE]  # Extract first nfs loadings

  # Varimax rotation
  rotated = varimax(loadings)
  A = unclass(rotated$loadings)

  # Eigenvalues = variance explained by each PC
  lambda = pc$sdev[1:nfs]^2

  # Scale loadings by sqrt of corresponding eigenvalues
  B = A / sqrt(matrix(rep(lambda, times = nrow(A)), ncol = ncol(A), byrow = TRUE))

  # Variance accounted by each component
  varP = lambda / sum(lambda) * 100  # Percentage of variance explained

  # Compute weighted contribution of each loading
  beta = abs(B %*% varP) / 100

  # Normalize weights
  w = beta / sum(beta)

  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)

  list(w = w[,1], s = s, lambda = lambda, B = B, beta = beta[,1])
}
#' @title PCA Weighting Method
