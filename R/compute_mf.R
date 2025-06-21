#' Compute fuzzy membership vector and return corresponding membership functions.
#'
#' \code{compute_mf} transforms a single indicator value into a fuzzy membership vector,
#' where each element represents the degree of membership to a specific evaluation level.
#' \code{compute_mf_funs} returns the list of membership functions for visualization purposes.
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
#' compute_mf(0.07, th)
#'
#' \dontrun{
#' mfs = compute_mf_funs(th)
#' plots = lapply(mfs, \(x) plot_mf(x, xlim = c(0, 0.6)))
#' gridExtra::grid.arrange(grobs = plots, nrow = 2)
#' }
#'
#'
#' @name compute_mf
NULL

#' @rdname compute_mf
#' @export
compute_mf_funs = function(thresholds) {
  # Generate membership functions based on thresholds (half trapezoids at both ends)
  # Input:
  #   thresholds: vector of intermediate thresholds (length >= 2)
  # Output:
  #   mfs: list of membership functions for all levels

  n = length(thresholds)
  if(n < 2) stop("thresholds must contain at least two values")

  mfs = list()

  # Right-half trapezoid for the first level
  mfs[[1]] = function(x) trap_mf(x, c(-Inf,-Inf,thresholds[1:2]))

  # Left-half trapezoid for the last level
  mfs[[n]] = function(x) trap_mf(x, c(thresholds[(n-1):n],Inf,Inf))

  if(n > 2) {
    mfs[2:(n-1)] = sapply(1:(n-2), function(i) {
      function(x) tri_mf(x, thresholds[i:(i + 2)])
    })
  }

  mfs
}


#' @rdname compute_mf
#' @export

compute_mf = function(x, thresholds) {
  # Convert a single indicator value into a fuzzy membership vector
  # Input:
  #   x: single indicator value
  #   thresholds: vector of intermediate thresholds (length >= 2)
  # Output:
  #   mv: membership vector

  n = length(thresholds)
  if(n < 2) stop("thresholds must contain at least two values")

  mfs = compute_mf_funs(thresholds)
  sapply(mfs, function(f) f(x))
}


