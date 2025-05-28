#' @title Combine Multiple Prediction Results
#'
#' @description Combines multiple prediction results (e.g., from grey prediction, time series, or machine learning models) into a single prediction using a similarity-based weighting approach, improving prediction accuracy.
#'
#' @param x Numeric vector, prediction results to be combined (length >= 2).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{a}: Numeric, the combined prediction value.
#'   \item \code{w}: Numeric vector, weights for each prediction in \code{x}, summing to 1.
#' }
#' @details
#' The function combines prediction results by constructing a similarity matrix based on cosine transformation of pairwise differences. Weights are derived from the principal eigenvector of the similarity matrix, ensuring predictions closer to each other have higher influence. For two predictions, equal weights (0.5, 0.5) are used. If all predictions are identical, equal weights are assigned.
#' Compatible with the \code{mathmodels} package for enhancing prediction models, including grey prediction, time series, or ensemble machine learning.
#' @examples
#' # Example: Combine three prediction results
#' preds = c(100, 102, 98)  # E.g., from grey prediction, ARIMA, or ML models
#' combine_preds(preds)

#' @export
combine_preds = function(x) {
  # Combine multiple prediction results
  Lx = length(x)
  if (Lx < 2) stop("Input x must have length >= 2.")
  if(Lx == 2) {    # For two predictions, use equal weights
    a = mean(x)
    w = c(1/2,1/2)
  }
  IndCom = combn(1:Lx,2)      # All pairwise combinations of element indices in x
  d = abs(x[IndCom[1,]] - x[IndCom[2,]])  # Calculate the distance between any two values
  maxd = max(d)
  idx = expand.grid(1:Lx,1:Lx)
  # Construct similarity matrix
  R = matrix(cos(pi*(x[idx[,1]]-x[idx[,2]]) / (2*maxd)),
             nrow = Lx, byrow = TRUE)
  rlt = eigen(R)
  w = rlt$vectors[,1] / sum(rlt$vectors[,1])
  a = sum(x * w)
  list(a = a, w = w)
}
