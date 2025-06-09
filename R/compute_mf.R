#' Compute fuzzy membership vector for a single indicator value.
#'
#' This function transforms a single indicator value into a fuzzy membership vector,
#' where each element represents the degree of membership to a specific evaluation level.
#' The first and last levels are modeled as half-trapezoidal functions (right and left side respectively),
#' while intermediate levels use triangular membership functions.
#'
#' @param x A numeric scalar representing the value of an indicator.
#' @param thresholds A numeric vector containing at least two threshold values that define the boundaries between evaluation levels.
#' @return A numeric vector with named elements indicating the membership degrees to each level.
#'
#' @examples
#' # Example: SO2 concentration = 0.07, thresholds = c(0.05, 0.15, 0.25, 0.5)
#' th = c(0.05, 0.15, 0.25, 0.5)
#' compute_mf(0.07, th)
#'
#' @export
compute_mf = function(x, thresholds) {
  n = length(thresholds)
  if(n < 2) stop("thresholds must contain at least two values")
  mv = numeric(n)

  # Right-half trapezoid for the first level
  mv[1] = (function(x) trap_mf(x, c(-Inf, -Inf, thresholds[1:2])))(x)

  # Left-half trapezoid for the last level
  mv[n] = (function(x) trap_mf(x, c(thresholds[(n-1):n], Inf, Inf)))(x)

  if(n == 2) return(mv)

  # Triangular membership for middle levels
  for(i in 1:(n - 2)) {
    mv[i + 1] = (function(x) tri_mf(x, thresholds[i:(i + 2)]))(x)
  }

  return(mv)
}
