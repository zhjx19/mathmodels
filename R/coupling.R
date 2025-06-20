#' Coupling and Coordination Analysis
#'
#' @description
#' Computes coupling degree, coordination index, and coupling coordination degree for multiple subsystems.
#'
#' @param u Data frame or matrix with normalized scores of subsystems as columns.
#' @param w Optional vector of weights for subsystems; defaults to equal weights if NULL.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{CD}{Coupling degree.}
#'   \item{CI}{Coordination index.}
#'   \item{CCD}{Coupling coordination degree.}
#' }
#'
#' @examples
#' # Sample normalized subsystem scores
#' df = data.frame(
#'   s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
#'   s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573)
#' )
#' # Equal weights
#' coupling(df)
#' # Custom weights
#' coupling(df, c(0.6, 0.4))
#'
#' @export
coupling = function(u, w = NULL) {
  # Computes coupling degree(CD), coordination index(CI), and coupling coordination degree(CCD)
  # u: data frame or matrix with normalized subsystem scores as columns
  # w: vector of subsystem weights
  p = ncol(u)
  if(is.null(w)) w = rep(1/p,p)
  row_prod = apply(u, 1, function(x) prod(pmax(x, 1e-10)))
  row_sum = pmax(rowSums(u), 1e-10)     # Prevent division by 0
  C = p * row_prod^(1/p) / row_sum
  TT = as.matrix(u) %*% w |> as.vector()
  D = sqrt(C * TT)
  data.frame(CD = C, CI = TT, CCD = D)
}
