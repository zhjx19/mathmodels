#' @title Defuzzification Methods for Fuzzy Comprehensive Evaluation
#'
#' @description Implements defuzzification methods for fuzzy evaluation vectors, including weighted average and maximum membership methods.
#'
#' @param mu Numeric vector, membership degrees for evaluation levels, in \eqn{[0, 1]}.
#' @param scores Numeric vector, scores corresponding to each evaluation level (e.g., c(100, 80, 60, 40) for "Excellent", "Good", "Fair", "Poor").
#' @param method Character, defuzzification method: "weighted_average", "max_membership", "centroid".
#' @return Numeric, defuzzified output value.
#'
#' @examples
#' # Example: Defuzzify fuzzy evaluation vectors for three schemes
#' mu = c(0.318, 0.351, 0.203, 0.128)
#' scores = c(30, 60, 75, 90)  # Scores for "Poor", "Fair", "Good", "Excellent"
#' defuzzify(mu, scores, method = "weighted_average")
#' defuzzify(mu, scores, method = "max_membership")
#' defuzzify(mu, scores, method = "centroid")

#' @export
defuzzify = function(mu, scores, method = "weighted_average") {
  switch(method,
         weighted_average = sum(mu * scores),
         max_membership = scores[mu == max(mu)],
         centroid = sum(mu * scores) / sum(mu)
  )
}
