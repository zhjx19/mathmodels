#' @title System Evaluation Functions for Coupling and Obstacle Analysis
#'
#' @description
#' These functions provide two key tools for system-level evaluation in multi-indicator systems:
#'
#' \itemize{
#'   \item \code{coupling_degree()}: Computes coupling degree, coordination index, and coupling coordination degree for subsystems.
#'   \item \code{obstacle_degree()}: Computes obstacle degree of each indicator to identify key constraints in the system.
#' }
#'
#' @param data A numeric matrix or data frame with normalized scores (usually in [0,1]) as columns.
#' @param w Optional vector of weights for indicators or subsystems; defaults to equal weights if NULL.
#'
#' @return
#' A list or data frame depending on the function:
#' \describe{
#'   \item{coupling_degree}{Data frame with columns:
#'     \itemize{
#'       \item \code{CD}: Coupling Degree (range 0-1)
#'       \item \code{CI}: Coordination Index (range 0-1)
#'       \item \code{CCD}: Coupling Coordination Degree (range 0-1)
#'     }
#'   }
#'   \item{obstacle_degree}{Data frame where each row sums to 100, showing percentage contribution of each indicator to total deviation.}
#' }
#'
#' @examples
#' # Sample normalized subsystem scores
#' df = data.frame(
#'   s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
#'   s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573)
#' )
#' # Coupling Degree Analysis
#' coupling_degree(df)        # Equal weights
#' coupling_degree(df, c(0.6, 0.4))
#' # Obstacle Degree Analysis
#' obstacle_degree(df)        # Equal weights
#' obstacle_degree(df, c(0.6, 0.4))
#'
#' @name system_evaluation
NULL

#' @rdname system_evaluation
#' @export
coupling_degree = function(data, w = NULL) {
  # Computes coupling degree(CD), coordination index(CI), and coupling coordination degree(CCD)
  # data: data frame or matrix with normalized subsystem scores as columns
  # w: vector of subsystem weights
  p = ncol(data)
  if(is.null(w)) w = rep(1/p,p)
  row_prod = apply(data, 1, function(x) prod(pmax(x, 1e-10)))
  row_sum = pmax(rowSums(data), 1e-10)   # Prevent division by 0
  C = p * row_prod^(1/p) / row_sum
  TT = as.matrix(data) %*% w |> as.vector()
  D = sqrt(C * TT)
  data.frame(CD = C, CI = TT, CCD = D)
}

#' @rdname system_evaluation
#' @export
obstacle_degree = function(data, w = NULL) {
  # Compute obstacle degree for each indicator
  # data: normalized data matrix or data frame
  # w: weights for indicators, default is equal weights
  n = nrow(data)
  m = ncol(data)
  if (is.null(w)) w = rep(1/m, m)
  diff_mat = 1 - as.matrix(data)
  weighted_diff = t(t(diff_mat) * w)
  res = 100 * weighted_diff / rowSums(weighted_diff)
  data.frame(res)
}
