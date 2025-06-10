#' Compute fuzzy membership vector and return corresponding membership functions.
#'
#' This function transforms a single indicator value into a fuzzy membership vector,
#' where each element represents the degree of membership to a specific evaluation level.
#' It also returns the list of membership functions for visualization purposes.
#'
#' @param x A numeric scalar representing the value of an indicator.
#' @param thresholds A numeric vector containing at least two threshold values that define
#'                   the boundaries between evaluation levels.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{mv}{A numeric vector, membership degrees to each level.}
#'     \item{mfs}{A list of functions, one per level, for plotting membership functions.}
#'   }
#'
#' @examples
#' # Example: SO2 concentration = 0.07, thresholds = c(0.05, 0.15, 0.25, 0.5)
#' th = c(0.05, 0.15, 0.25, 0.5)
#' res = compute_mf(0.07, th)
#' res$mv
#'
#' \dontrun{
#' plots = lapply(res$mfs, \(x) plot_mf(x, xlim = c(0, 0.6)))
#' gridExtra::grid.arrange(grobs = plots, nrow = 2)
#' }
#'
#' @export
compute_mf = function(x, thresholds) {
  # Based on triangular membership (with half-trapezoidal at both ends),
  # convert a single indicator value into a membership vector over evaluation levels.
  # Input:
  #   x: single indicator value
  #   thresholds: vector of intermediate thresholds (length >= 2)
  # Output:
  #   mv: membership vector
  #   mfs: list of membership functions for each level

  n = length(thresholds)
  if(n < 2) stop("thresholds must contain at least two values")
  mv = numeric(n)
  mfs = list()

  # Right-half trapezoid for the first level
  mfs[[1]] = function(x) trap_mf(x, c(-Inf,-Inf,thresholds[1:2]))
  mv[1] = (mfs[[1]])(x)

  # Left-half trapezoid for the last level
  mfs[[n]] = function(x) trap_mf(x, c(thresholds[(n-1):n],Inf,Inf))
  mv[n] = (mfs[[n]])(x)

  if(n == 2) return(list(mv = mv, mfs = mfs))

  # Triangular membership for middle levels
  for(i in 1:(n-2)) {
    mfs[[i+1]] = function(x) tri_mf(x, thresholds[i:(i+2)])
    mv[i+1] = (mfs[[i+1]])(x)
  }

  return(list(mv = mv, mfs = mfs))
}
