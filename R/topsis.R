#' TOPSIS Method for Multi-Criteria Decision Making
#'
#' @description
#' Implements the Technique for Order of Preference by Similarity to Ideal Solution
#' (TOPSIS) method to rank alternatives based on multiple criteria. The decision matrix
#' \code{A} must contain positive indicators (larger values are better) or be preprocessed
#' using functions like \code{\link{rescale}}, \code{\link{normalize}},
#' \code{\link{rescale_middle}}, or \code{\link{rescale_interval}} from the
#' \pkg{mathmodels} package. The function normalizes the matrix, applies weights, and
#' computes relative closeness to the ideal solution.
#'
#' @param A Numeric matrix, the decision matrix where rows are alternatives and columns
#' are criteria. All criteria must be positive indicators (larger is better) or preprocessed
#' to satisfy this condition.
#' @param w Numeric vector, weights for each criterion. Must be non-negative and sum to
#' a positive value.
#'
#' @return A named numeric vector of the same length as the number of rows in \code{A},
#' containing relative closeness scores in [0, 1]. Higher values indicate better alternatives.
#' Names are taken from \code{rownames(A)} or default to "Sample1", "Sample2", etc.
#'
#' @details
#' The TOPSIS method ranks alternatives by:
#' \enumerate{
#'   \item Normalizing the decision matrix using L2 norm (similar to \code{\link{normalize}}).
#'   \item Applying weights to form a weighted normalized matrix.
#'   \item Identifying positive and negative ideal solutions (column max/min).
#'   \item Computing Euclidean distances to ideal solutions.
#'   \item Calculating relative closeness as \code{S0 / (S0 + Sstar)}, where \code{S0}
#'   is the distance to the negative ideal and \code{Sstar} is the distance to the positive ideal.
#' }
#' Since \code{A} must contain positive indicators, use \code{\link{rescale}} for Min-Max
#' normalization, \code{\link{rescale_middle}} for centered indicators, or
#' \code{\link{rescale_interval}} for interval indicators before calling \code{topsis()}.
#' NA values in \code{A} may affect results and trigger a warning.
#'
#' @references
#' Hwang, C. L., & Yoon, K. (1981). Multiple Attribute Decision Making: Methods and Applications.
#' The \code{mathmodels_penguins} dataset is derived from the \pkg{palmerpenguins} package.
#'
#' @export
#' @examples
#' A = data.frame(
#'   X1 = c(2, 5, 3),  # "+"
#'   X2 = c(8, 1, 6)   # "+"
#' )
#' w = c(0.6, 0.4)
#' topsis(A, w)

topsis = function(A, w) {
  # Implements TOPSIS method, A is decision matrix, must be positive indicators or preprocessed
  # to satisfy this condition. w is weight vector
  n = nrow(A)
  B = apply(A, 2, normalize)         # Normalize the decision matrix
  C = t(t(B) * w)                    # Weighted normalized matrix, or
  # C = B * matrix(rep(w, each = n), nrow = n)
  Cstar = apply(C, 2, max)           # Column max for positive ideal solution
  C0 = apply(C, 2, min)              # Column min for negative ideal solution
  # Distance to positive ideal solution
  Sstar = apply(C, 1, \(x) norm(x - Cstar, "2"))
  # Distance to negative ideal solution
  S0 = apply(C, 1, \(x) norm(x - C0, "2"))
  S0 / (S0 + Sstar)                  # Relative closeness
}
