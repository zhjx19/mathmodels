#' Grey Relational Analysis Functions
#'
#' @description
#' A collection of functions for performing grey relational analysis, including
#' calculation of grey correlation degree and evaluation based on grey correlation.
#' These functions are designed for decision-making and data analysis by measuring
#' the relational degree between sequences.
#'
#' @param ref Numeric vector, the reference sequence for `grey_corr`.
#' @param cmp Numeric matrix or data frame, the comparison sequences for `grey_corr`.
#' @param rho Numeric scalar, the distinguishing coefficient (default = 0.5).
#' @param w Numeric vector, weights for weighted correlation (default = equal weights).
#' @param X Numeric matrix or data frame, the decision matrix for `grey_corr_topsis`.
#' @param index Character vector indicating indicator direction:
#'   Use `"+"` for positive indicators (higher is better),
#'   `"-"` for negative indicators (lower is better),
#'   and `NA` for already rescaled indicators (no rescaling will be applied).
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
#'   \item{grey_corr}{Computes the grey correlation degree between a reference sequence (`ref`)
#'                    and comparison sequences (`cmp`) using the distinguishing coefficient (`rho`)
#'                    and optional weights (`w`).}
#'   \item{grey_corr_topsis}{Evaluates a decision matrix (`X`) by normalizing it,
#'                           applying weights (`w`), computing grey correlation with the ideal sequence.
#'                           Direction of indicators can be specified via `index`.}
#' }
#'
#' @examples
#' # Grey correlation degree
#' ref = c(0.9, 0.8, 0.7)
#' cmp = data.frame(
#'   x1 = c(0.9, 0.7, 0.8),
#'   x2 = c(0.8, 0.9, 0.7),
#'   x3 = c(0.7, 0.8, 0.9)
#' )
#' grey_corr(ref, cmp, rho = 0.5)
#'
#' # Grey correlation evaluation
#' X = data.frame(x1 = c(8, 7, 6), x2 = c(150, 180, 200), x3 = c(60, 80, 100))
#' w = c(0.3, 0.4, 0.3)
#' idx = c("+", "+", "+")
#' grey_corr_topsis(X, w, idx, rho = 0.5)
#'
#' @name grey_analysis
NULL

#' @rdname grey_analysis
#' @export
grey_corr = function(ref, cmp, rho = 0.5, w = NULL) {
  # Compute grey correlation degree
  # ref: reference sequence, cmp: comparison sequences, rho: distinguishing coefficient
  # w: weights for weighted correlation

  if(!is.numeric(ref)) stop("ref must be a numeric vector.")
  if(!is.matrix(cmp) && !is.data.frame(cmp))
    stop("cmp must be a matrix or data frame.")
  if(length(ref) != nrow(cmp))
    stop("Length of ref must equal number of rows in cmp.")
  if(!is.null(w) && length(w) != nrow(cmp))
    stop("w must have length equal to nrow(cmp).")

  p = nrow(cmp)
  if(is.null(w)) w = rep(1/p, p)
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

  if(!is.matrix(X) && !is.data.frame(X))
    stop("X must be a matrix or data frame.")
  m = ncol(X)
  if(is.null(index)) index = rep(NA, m)
  if(length(index) != m)
    stop("index must have length equal to ncol(X).")
  if(is.null(w)) w = rep(1/m, m)
  if(length(w) != m)
    stop("w must have length equal to ncol(X).")
  if(any(w < 0) || sum(w) == 0)
    stop("w must be non-negative and sum to a positive value.")

  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], function(x) rescale(x, type = "-"))
  # Weighted normalized matrix
  C = as.matrix(X) %*% diag(w)
  Cstar = apply(C, 2, max)
  C0 = apply(C, 2, min)
  Sstar = grey_corr(Cstar, C, rho)
  S0 = grey_corr(C0, C, rho)
  Sstar / (S0 + Sstar)
}
