# Evaluation Algorithms — Fuzzy Evaluation
#
# Includes: membership functions, fuzzy_eval, compute_mf, and defuzzify
#


# ============================================================================
# 隶属函数 (Membership Functions)
# ============================================================================

#' @title Membership Functions for Fuzzy Logic
#'
#' @description
#' A collection of functions to compute membership values for various fuzzy sets, including triangular, trapezoidal, Gaussian, generalized bell, two-parameter Gaussian, sigmoid, difference of sigmoids, product of sigmoids, Z-shaped, PI-shaped, and S-shaped membership functions. Includes a function to visualize membership functions using ggplot2. These are designed for evaluation models in mathematical modeling, compatible with `fuzzy_eval` in the `mathmodels` package.
#'
#' @param x Numeric vector, input values for which to compute membership.
#' @param params Numeric vector, parameters defining the membership function:
#' \itemize{
#'   \item For `tri_mf`: `c(a, b, c)`, where `a <= b <= c` (left base, peak, right base).
#'   \item For `trap_mf`: `c(a, b, c, d)`, where `a <= b <= c <= d` (left base, left top, right top, right base).
#'   \item For `gauss_mf`: `c(sigma, c)`, where `sigma > 0` (spread, center).
#'   \item For `gbell_mf`: `c(a, b, c)`, where `a > 0`, `b > 0` (width, shape, center).
#'   \item For `gauss2mf`: `c(s1, c1, s2, c2)`, where `s1 > 0`, `s2 > 0` (left spread, left center, right spread, right center).
#'   \item For `sigmoid_mf`: `c(a, b)`, where `a > 0` (slope, inflection point).
#'   \item For `dsigmoid_mf`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0` (slopes and inflection points for two sigmoids).
#'   \item For `psigmoid_mf`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0` (slopes and inflection points for two sigmoids).
#'   \item For `z_mf`: `c(a, b)`, where `a < b` (left base, right base).
#'   \item For `pi_mf`: `c(a, b, c, d)`, where `a < b < c < d` (left base, left shoulder, right shoulder, right base).
#'   \item For `s_mf`: `c(a, b)`, where `a < b` (left base, right base).
#' }
#' @param mf Function, a membership function with fixed parameters (e.g., `function(x) tri_mf(x, c(2, 5, 8))`).
#' @param xlim Numeric vector of length 2, x-axis limits for plotting (default `c(0, 10)`).
#' @param main Character, plot title (default `NULL`, no title).
#'
#' @return
#' \itemize{
#'   \item For membership functions (`tri_mf`, `trap_mf`, `gauss_mf`, `gbell_mf`, `gauss2mf`, `sigmoid_mf`, `dsigmoid_mf`, `psigmoid_mf`, `z_mf`, `pi_mf`, `s_mf`): A numeric vector of membership values in \code{[0, 1]}, same length as `x`.
#'   \item For `plot_mf`: A ggplot2 object, plotting the membership function.
#' }
#'
#' @details
#' These functions support evaluation models in mathematical modeling:
#' \itemize{
#'   \item `tri_mf`: Triangular membership, linear rise from `a` to `b` (peak) and fall to `c`.
#'   \item `trap_mf`: Trapezoidal membership, linear rise from `a` to `b`, plateau from `b` to `c`, fall to `d`.
#'   \item `gauss_mf`: Gaussian membership, bell-shaped curve centered at `c` with spread `sigma`.
#'   \item `gbell_mf`: Generalized bell membership, bell-shaped curve with width `a`, shape `b`, and center `c`.
#'   \item `gauss2mf`: Two-parameter Gaussian membership, combining two Gaussians with spreads `s1`, `s2` and centers `c1`, `c2`.
#'   \item `sigmoid_mf`: Sigmoid membership, S-shaped curve with slope `a` and inflection point `b`.
#'   \item `dsigmoid_mf`: Difference of two sigmoids, combining slopes `a1`, `a2` and inflection points `c1`, `c2`.
#'   \item `psigmoid_mf`: Product of two sigmoids, combining slopes `a1`, `a2` and inflection points `c1`, `c2`.
#'   \item `z_mf`: Z-shaped membership, decreasing from 1 at `a` to 0 at `b`.
#'   \item `pi_mf`: PI-shaped membership, rising from `a` to `b`, plateau from `b` to `c`, falling to `d`.
#'   \item `s_mf`: S-shaped membership, increasing from 0 at `a` to 1 at `b`.
#'   \item `plot_mf`: Plots a membership function over `xlim` using ggplot2, suitable for tidyverse workflows.
#' }
#' Membership values can be used to construct fuzzy evaluation matrices for `fuzzy_eval`.
#' Implemented in base R, except `plot_mf`, which requires ggplot2.
#'
#' @examples
#' # Define input values
#' x = 0:10
#'
#' # Triangular membership
#' tri_mf(x, params = c(3, 6, 8))
#'
#' # Trapezoidal membership
#' trap_mf(x, params = c(1, 5, 7, 8))
#'
#' # Gaussian membership
#' gauss_mf(x, params = c(2, 5))
#'
#' # Generalized bell membership
#' gbell_mf(x, params = c(2, 4, 6))
#'
#' # Two-parameter Gaussian membership
#' gauss2mf(x, params = c(1, 3, 3, 4))
#'
#' # Sigmoid membership
#' sigmoid_mf(x, params = c(2, 4))
#'
#' # Difference of sigmoids membership
#' dsigmoid_y = dsigmoid_mf(x, params = c(5, 2, 5, 7))
#'
#' # Product of sigmoids membership
#' psigmoid_mf(x, params = c(2, 3, -5, 8))
#'
#' # Z-shaped membership
#' z_mf(x, params = c(3, 7))
#'
#' # PI-shaped membership
#' pi_mf(x, params = c(1, 4, 5, 10))
#'
#' # S-shaped membership
#' s_mf(x, params = c(1, 8))
#'
#' \dontrun{
#' # Visualize membership functions
#' plot_mf(\(x) tri_mf(x, c(3, 6, 8)), main = "Triangular MF")
#' plot_mf(\(x) trap_mf(x, c(1, 5, 7, 8)), main = "Trapezoidal MF")
#' plot_mf(\(x) gauss_mf(x, c(2, 5)), main = "Gaussian MF")
#' plot_mf(\(x) gbell_mf(x, c(2, 4, 6)), main = "Generalized Bell MF")
#' plot_mf(\(x) gauss2mf(x, c(1, 3, 3, 4)), main = "Two-Parameter Gaussian MF")
#' plot_mf(\(x) sigmoid_mf(x, c(2, 4)), main = "Sigmoid MF")
#' plot_mf(\(x) dsigmoid_mf(x, c(5, 2, 5, 7)), main = "Difference of Sigmoids MF")
#' plot_mf(\(x) psigmoid_mf(x, c(2, 3, -5, 8)), main = "Product of Sigmoids MF")
#' plot_mf(\(x) z_mf(x, c(3, 7)), main = "Z-Shaped MF")
#' plot_mf(\(x) pi_mf(x, c(1, 4, 5, 10)), main = "PI-Shaped MF")
#' plot_mf(\(x) s_mf(x, c(1, 8)), main = "S-Shaped MF")
#' }
#'
#' @name membership
NULL

#' @rdname membership
#' @export
tri_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 3)
    stop("params must be a numeric vector of length 3.")
  a = params[1]
  b = params[2]
  C = params[3]
  y = rep(0, length(x))
  y[x >= a & x <= b] = (x[x >= a & x <= b] - a) / (b - a)
  y[x > b & x <= C] = (C - x[x > b & x <= C]) / (C - b)
  y
}

#' @rdname membership
#' @export
trap_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  y = rep(0, length(x))
  y[x >= a & x <= b] = (x[x >= a & x <= b] - a) / (b - a)
  y[x > b & x < C] = 1
  y[x >= C & x <= d] = (d - x[x >= C & x <= d]) / (d - C)
  y
}

#' @rdname membership
#' @export
gauss_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  sigma = params[1]
  C = params[2]
  exp(-(x - C)^2 / (2 * sigma^2))
}

#' @rdname membership
#' @export
gbell_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 3)
    stop("params must be a numeric vector of length 3.")
  a = params[1]
  b = params[2]
  C = params[3]
  1 / (1 + abs((x - C) / a)^(2*b))
}

#' @rdname membership
#' @export
gauss2mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  s1 = params[1]
  c1 = params[2]
  s2 = params[3]
  c2 = params[4]
  y = rep(1, length(x))
  y[x < c1] = exp(-(x[x<c1] - c1)^2 / (2 * s1^2))
  y[x > c2] = exp(-(x[x>c2] - c2)^2 / (2 * s2^2))
  y
}

# Sigmoid helper (internal)
sigmoid = function(x, a, b) 1 / (1 + exp(-a * (x - b)))

#' @rdname membership
#' @export
sigmoid_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a1 = params[1]
  c1 = params[2]
  sigmoid(x, a1, c1)
}

#' @rdname membership
#' @export
dsigmoid_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a1 = params[1]
  c1 = params[2]
  a2 = params[3]
  c2 = params[4]
  sigmoid(x, a1, c1) - sigmoid(x, a2, c2)
}

#' @rdname membership
#' @export
psigmoid_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a1 = params[1]
  c1 = params[2]
  a2 = params[3]
  c2 = params[4]
  sigmoid(x, a1, c1) * sigmoid(x, a2, c2)
}

#' @rdname membership
#' @export
z_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  y = rep(0, length(x))
  y[x <= a] = 1
  cond = x > a & x <= (a+b)/2
  y[cond] = 1 - 2 * ((x[cond] - a) / (b - a))^2
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 2 * ((x[cond2] - b) / (b - a))^2
  y
}

#' @rdname membership
#' @export
pi_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  y = rep(0, length(x))
  y[x >= b & x <= C] = 1
  cond = x > a & x <= (a+b)/2
  y[cond] = 2 * ((x[cond] - a) / (b - a))^2
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 1 - 2 * ((x[cond2] - b) / (b - a))^2
  cond3 = x > C & x <= (C+d)/2
  y[cond3] = 1 - 2 * ((x[cond3] - C) / (d - C))^2
  cond4 = x > (C+d)/2 & x < d
  y[cond4] = 2 * ((x[cond4] - d) / (d - C))^2
  y
}

#' @rdname membership
#' @export
s_mf = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  y = rep(0, length(x))
  y[x >= b] = 1
  cond = x > a & x <= (a+b)/2
  y[cond] = 2 * ((x[cond] - a) / (b - a))^2
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 1 - 2 * ((x[cond2] - b) / (b - a))^2
  y
}

#' @rdname membership
#' @export
#' @importFrom rlang .data
plot_mf = function(mf, xlim = c(0, 10), main = NULL) {
  if(!is.function(mf)) stop("mf must be a function.")
  data = data.frame(x = seq(xlim[1], xlim[2], length.out = 200),
                     y = mf(seq(xlim[1], xlim[2], length.out = 200)))
  p = ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1.2) +
    ggplot2::ylim(0, 1.05) + ggplot2::theme_minimal()
  if(!is.null(main)) p = p + ggplot2::labs(title = main)
  p
}


# ============================================================================
# compute_mf — 指标值 → 模糊隶属向量
# ============================================================================

#' Compute fuzzy membership vector and return corresponding membership functions.
#'
#' `compute_mf` transforms a single indicator value into a fuzzy membership vector,
#' where each element represents the degree of membership to a specific evaluation level.
#' `compute_mf_funs` returns the list of membership functions for visualization purposes.
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
#' @name compute_mf
NULL

#' @rdname compute_mf
#' @export
compute_mf_funs = function(thresholds) {
  if(!is.numeric(thresholds) || is.matrix(thresholds))
    stop("thresholds must be a numeric vector.")
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
  if(!is.numeric(x) || length(x) != 1)
    stop("x must be a numeric scalar.")
  if(!is.numeric(thresholds) || is.matrix(thresholds))
    stop("thresholds must be a numeric vector.")
  n = length(thresholds)
  if(n < 2) stop("thresholds must contain at least two values")

  mfs = compute_mf_funs(thresholds)
  sapply(mfs, function(f) f(x))
}


# ============================================================================
# fuzzy_eval — 模糊综合评价
# ============================================================================

#' Fuzzy Comprehensive Evaluation
#'
#' @description
#' Performs fuzzy comprehensive evaluation using different fuzzy composition operators
#' to combine factor weights with a fuzzy evaluation matrix. Suitable for multi-criteria
#' decision analysis with weights from methods like AHP, entropy, CRITIC, CV, or PCA.
#'
#' @param w Numeric vector, factor weights (e.g., from `combine_weights`).
#' @param R Numeric matrix, fuzzy evaluation matrix with columns as factors and rows
#' as evaluation grades. Values should be in \code{[0, 1]}.
#' @param type Integer or character (1-5), specifying the fuzzy composition operator:
#' \itemize{
#'   \item 1: Min-max (main factor decisive).
#'   \item 2: Product-max (main factor prominent).
#'   \item 3: Weighted sum (additive average).
#'   \item 4: Bounded sum of mins (min-sum bounded).
#'   \item 5: Normalized min-sum (balanced average).
#' }
#'
#' @return A numeric vector of normalized comprehensive evaluation results, summing to 1.
#'
#' @details
#' The function computes a fuzzy comprehensive evaluation vector `B` based on
#' the weight vector `w` and fuzzy evaluation matrix `R`. Five composition
#' operators are supported:
#' \itemize{
#'   \item Type 1 (min-max): `max(min(w, R[j,]))`, emphasizes the main factor.
#'   \item Type 2 (product-max): `max(w * R[j,])`, highlights the main factor.
#'   \item Type 3 (weighted sum): `sum(w * R[j,])`, additive average.
#'   \item Type 4 (bounded sum): `min(1, sum(min(w, R[j,])))`, bounds the sum of mins.
#'   \item Type 5 (normalized min-sum): `sum(min(w, R[j,]/sum(R[j,])))`, balanced average.
#' }
#' The output `B` is normalized to sum to 1. If the sum is zero, an error is thrown.
#' Uses base R for lightweight implementation.
#'
#' @examples
#' w = c(0.3, 0.3, 0.3, 0.1)  # weights (e.g., from AHP or entropy)
#'
#' # fuzzy evaluation matrix (3 grades for 4 factors)
#' R = matrix(c(0.8, 0.7, 0.6, 0.7,
#'              0.1, 0.2, 0.2, 0.1,
#'              0.1, 0.1, 0.2, 0.2), nrow = 3, byrow = TRUE)
#' # Apply fuzzy comprehensive evaluation
#' fuzzy_eval(w, R, type = 3)  # Weighted sum
#'
#' @export
fuzzy_eval = function(w, R, type) {
  if(!is.numeric(w) || is.matrix(w))
    stop("w must be a numeric vector.")
  if(!is.matrix(R))
    stop("R must be a matrix.")
  if(length(w) != ncol(R))
    stop("Length of w must equal ncol(R).")
  type = as.character(type)
  if(!type %in% as.character(1:5))
    stop("type must be 1, 2, 3, 4, or 5.")

  m = nrow(R)
  n = ncol(R)
  B = rep(0, m)
  for(j in 1:m) {
    switch(type,
           # Min-max, main factor decisive
           "1" = {
             B[j] = max(apply(rbind(w, R[j,]), 2, min))
           },
           #  Product-max, main factor prominent
           "2" = {
             B[j] = max(w * R[j,])
           },
           # Weighted sum, additive average
           "3" = {
             B[j] = sum(w * R[j,])
           },
           # Bounded sum of mins
           "4" = {
             B[j] = min(c(1, sum(apply(rbind(w, R[j,]), 2, min))))
           },
           # Normalized min-sum, balanced average
           "5" = {
             r0 = sum(R[j,])
             B[j] = sum(apply(rbind(w, R[j,] / r0), 2, min))
           }
    )
  }
  if(sum(B) == 0) stop("Normalization failed: sum of B is zero.")
  B / sum(B)      # Normalize
}


# ============================================================================
# defuzzify — 去模糊化
# ============================================================================

#' @title Defuzzification Methods for Fuzzy Comprehensive Evaluation
#'
#' @description Implements defuzzification methods for fuzzy evaluation vectors, including weighted average and maximum membership methods.
#'
#' @param mu Numeric vector, membership degrees for evaluation levels, in \code{[0, 1]}.
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
  if(!is.numeric(mu) || is.matrix(mu))
    stop("mu must be a numeric vector.")
  if(!is.numeric(scores) || is.matrix(scores))
    stop("scores must be a numeric vector.")
  if(length(mu) != length(scores))
    stop("mu and scores must have the same length.")
  switch(method,
         weighted_average = sum(mu * scores),
         max_membership = scores[mu == max(mu)],
         centroid = sum(mu * scores) / sum(mu)
  )
}
