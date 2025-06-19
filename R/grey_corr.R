#' Grey Relational Analysis Functions
#'
#' @description
#' A collection of functions for performing grey relational analysis, including
#' calculation of grey correlation degree and evaluation based on grey correlation.
#' These functions are designed for decision-making and data analysis by measuring
#' the relational degree between sequences.
#'
#' @param ref Numeric vector, the reference sequence for \code{grey_corr}.
#' @param cmp Numeric matrix or data.frame, the comparison sequences for \code{grey_corr}.
#' @param rho Numeric scalar, the distinguishing coefficient (default = 0.5).
#' @param w Numeric vector, weights for weighted correlation (default = equal weights).
#' @param A Numeric matrix or data.frame, the decision matrix for \code{grey_corr_eval}.
#'
#' @return
#' \itemize{
#'   \item \code{grey_corr}: Numeric vector, the grey correlation degree for each comparison sequence.
#'   \item \code{grey_corr_eval}: Numeric vector, normalized evaluation scores in [0, 100].
#' }
#'
#' @details
#' These functions implement grey relational analysis for evaluating relationships
#' between sequences or decision alternatives:
#' \itemize{
#'   \item \code{grey_corr}: Computes the grey correlation degree between a reference
#'   sequence (\code{ref}) and comparison sequences (\code{cmp}) using the distinguishing
#'   coefficient (\code{rho}) and optional weights (\code{w}).
#'   \item \code{grey_corr_eval}: Evaluates a decision matrix (\code{A}) by normalizing
#'   it, applying weights (\code{w}), computing grey correlation with the ideal sequence,
#'   and scaling results to [0, 100].
#' }
#'
#' @examples
#' # Grey correlation degree
#' ref = 1:3
#' cmp = data.frame(x1 = c(1, 2, 4), x2 = c(2, 3, 5))
#' grey_corr(ref, cmp, rho = 0.5)
#'
#' # Grey correlation evaluation#'
#' w = c(0.4, 0.6)
#' grey_corr_eval(cmp, w, rho = 0.5)
#'
#' @name grey_analysis
NULL

#' @rdname grey_analysis
#' @export
grey_corr = function(ref, cmp, rho = 0.5, w = NULL) {
  # Compute grey correlation degree
  # ref: reference sequence, cmp: comparison sequences, rho: distinguishing coefficient
  # w: weights for weighted correlation
  n = nrow(cmp)
  if(is.null(w)) w = rep(1/n, n)
  t = apply(cmp, 2, \(x) x - ref)
  min2 = min(apply(abs(t), 2, min))       # Minimum difference
  max2 = max(apply(abs(t), 2, max))       # Maximum difference
  eta = (min2 + rho*max2) / (abs(t) + rho*max2)  # Correlation coefficient
  r = w %*% eta
  r[1,]
}

#' @rdname grey_analysis
#' @export
grey_corr_eval = function(A, w, rho = 0.5) {
  # Perform grey correlation evaluation
  # A: decision matrix, w: weights for indicators
  # rho: distinguishing coefficient
  B = apply(A, 2, \(x) x / norm(x, "2"))  # Normalize matrix
  C = apply(B, 1, \(x) w * x) |> t()      # Weighted normalized matrix
  ref = apply(C, 2, max)
  grey_corr(ref, t(C), rho)
}
