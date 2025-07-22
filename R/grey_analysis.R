#' Grey Relational Analysis Functions
#'
#' @description
#' A collection of functions for performing grey relational analysis, including
#' calculation of grey correlation degree and evaluation based on grey correlation.
#' These functions are designed for decision-making and data analysis by measuring
#' the relational degree between sequences.
#'
#' @param ref Numeric vector, the reference sequence for \code{grey_corr}.
#' @param cmp Numeric matrix or data frame, the comparison sequences for \code{grey_corr}.
#' @param rho Numeric scalar, the distinguishing coefficient (default = 0.5).
#' @param w Numeric vector, weights for weighted correlation (default = equal weights).
#' @param X Numeric matrix or data frame, the decision matrix for \code{grey_corr_topsis}.
#' @param index Character vector indicating indicator direction:
#'   Use \code{"+"} for positive indicators (higher is better),
#'   \code{"-"} for negative indicators (lower is better),
#'   and \code{NA} for already rescaled indicators (no rescaling will be applied).
#'   If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no rescaling is performed.
#'
#' @return
#' \describe{
#'   \item{grey_corr}{Returns a numeric vector of grey correlation degrees for each comparison sequence.}
#'   \item{grey_corr_topsis}{Returns a numeric vector of relative closeness (grey correlation degrees).}
#' }
#'
#' @details
#' These functions implement grey relational analysis for evaluating relationships
#' between sequences or decision alternatives:
#' \describe{
#'   \item{grey_corr}{Computes the grey correlation degree between a reference sequence (\code{ref})
#'                    and comparison sequences (\code{cmp}) using the distinguishing coefficient (\code{rho})
#'                    and optional weights (\code{w}).}
#'   \item{grey_corr_topsis}{Evaluates a decision matrix (\code{X}) by normalizing it,
#'                           applying weights (\code{w}), computing grey correlation with the ideal sequence.
#'                           Direction of indicators can be specified via \code{index}.}
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
#' idx = c("+", "-")
#' grey_corr_topsis(cmp, w, idx, rho = 0.5)
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
grey_corr_topsis = function(X, w, index = NULL, rho = 0.5) {
  # Perform grey correlation TOPSIS evaluation
  # X: decision matrix
  # w: weights for indicators, defaults to equal weights
  # index: direction of each indicator, "+" for positive, "-" for negative, NA means no rescaling
  # rho: distinguishing coefficient
  m = ncol(X)
  if(is.null(index)) index = rep(NA, m)
  if(is.null(w)) w = rep(1/m, m)
  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], function(x) rescale(x, type = "-"))
  # Weighted normalized matrix
  C = as.matrix(X) %*% diag(w)
  Cstar = apply(C, 2, max)
  C0 = apply(C, 2, min)
  Sstar = grey_corr(Cstar, t(C), rho)
  S0 = grey_corr(C0, t(C), rho)
  Sstar / (S0 + Sstar)
}
