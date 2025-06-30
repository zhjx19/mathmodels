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
#' @param varvarimax Whether to perform Varimax rotation, default is TRUE.
#'
#' @param method Weighting Method, "abs" (default, |a_{ji}|) or "squared" (a_{ji}^2)
#'
#' @return A list containing:
#' \item{w}{Numeric vector of normalized weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample.}
#'
#' \item{lambda}{Eigenvalues of principal components (explained variance).}
#'
#' \item{varP}{Proportion of variance explained by selected PCs.}
#'
#' @export
#' @examples
#' # Example: Using PCA to compute indicator weights
#' ind = c("+","+","-","-")
#' pca_weight(iris[1:10, 1:4], ind, nfs = 2)

pca_weight = function(X, index = NULL, nfs = NULL, varimax = TRUE,
                      method = "abs") {
  # Compute indicator weights using Principal Component Analysis (PCA)
  # X: original data matrix, rows are samples, columns are indicators
  # index: direction of each indicator, "+" for positive, "-" for negative, `NA` for already normalized indicators
  # nfs: numbers of PCs to extract
  # varimax: whether to perform Varimax rotation
  # method: Weighting Method, "abs" (default, |a_{ji}|) or "squared" (a_{ji}^2)
  # w: returned normalized weights, lambda: eigenvalues
  # s: scores for each sample

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Standardize data
  idx = !is.na(index)
  X[idx] = scale(X[idx])
  X[,neg] = -X[,neg]

  # Perform PCA with all components and varimax rotation
  if(is.null(nfs)) nfs = ncol(X)

  # Replace psych::principal with prcomp + varimax
  pc = prcomp(X, center = FALSE, scale. = FALSE) # Already scaled above
  loadings = pc$rotation[, 1:nfs, drop = FALSE]  # Extract first nfs loadings

  # Varimax rotation
  if(varimax && nfs > 1) {
    rotated = varimax(loadings)
    A = unclass(rotated$loadings)
  } else {
    A = loadings
  }

  # Eigenvalues = variance explained by each PC
  eigenvalues = pc$sdev ^ 2
  lambda = eigenvalues[1:nfs]

  # Variance accounted by each component (proportion for selected components)
  varP = lambda / sum(eigenvalues)  # Proportion of variance explained by selected PCs

  # Compute weighted contribution of each loading
  beta = switch(method,
                "abs" = abs(A) %*% varP,
                "squared" = (A^2) %*% varP)
  # Normalize weights
  w = beta / sum(beta)

  # Compute sample scores
  s = as.vector(as.matrix(X) %*% w)  # Use standardized data

  list(w = as.vector(w), s = s, cum_contrib = sum(varP), loading = A,
       lambda = lambda, varP = varP)
}
#' @title PCA Weighting Method
