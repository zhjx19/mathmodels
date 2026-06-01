#' Coefficient of Variation Weighting
#'
#' @description
#' Computes weights for indicators using the Coefficient of Variation (CV) method.
#' Weights are derived by normalizing the CV (standard deviation divided by mean)
#' for each indicator.
#'
#' @param X Numeric matrix or data frame with positive indicator data.
#'
#' @return Numeric vector of weights for the indicators, summing to 1.
#'
#' @details
#' The `cv_weight` function calculates weights using the CV method.
#' For each column in `data`, the CV is computed as the standard deviation
#' divided by the mean. Weights are obtained by normalizing the CVs to sum to 1.
#' This lightweight implementation uses base R and assumes all columns are numeric
#' indicators.
#'
#' @examples
#' X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
#' cv_weight(X)
#'
#' @export
cv_weight = function(X) {
  # Compute weights using Coefficient of Variation
  # X: positive indicator data

  if(!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(any(apply(X, 2, function(x) any(x <= 0, na.rm = TRUE))))
    stop("All values in X must be positive for CV weighting.")

  # Calculate mean and standard deviation
  means = apply(X, 2, mean, na.rm = TRUE)
  sds = apply(X, 2, sd, na.rm = TRUE)
  cv = sds / means     # Compute coefficient of variation
  # Normalize to weights
  cv / sum(cv)
}
