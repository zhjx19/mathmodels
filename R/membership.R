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
  # Compute triangular membership function
  # x: input values, params: c(a, b, c) where a <= b <= c
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
  # Compute trapezoidal membership function
  # x: input values, params: c(a, b, c, d) where a <= b <= c <= d
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
  # Compute Gaussian membership function
  # x: input values, params: c(sigma, c) where sigma > 0
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
  # Compute generalized bell membership function
  # x: input values, params: c(a, b, c) where a > 0, b > 0, c > 0
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
  # Compute Gaussian membership function with two parameters
  # x: input values, params: c(s1, c1, s2, c2) where s1 > 0, s2 > 0
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

sigmoid = function(x, a, b) 1 / (1 + exp(-a * (x - b)))

#' @rdname membership
#' @export
sigmoid_mf = function(x, params) {
  # Compute sigmoid membership function
  # x: input values, params: c(a, b) where a > 0
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
  # Compute difference of sigmoid membership function
  # x: input values, params: c(a1, c1, a2, c2) where a1 > 0, a2 > 0
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
  # Compute product of sigmoid membership function
  # x: input values, params: c(a1, c1, a2, c2) where a1 > 0, a2 > 0
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
  # Compute Z-shaped membership function
  # x: input values, params: c(a, b) where a < b
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  y = rep(0, length(x))
  y[x <= a] = 1
  # Left half: a < x <= (a+b)/2
  cond = x > a & x <= (a+b)/2
  y[cond] = 1 - 2 * ((x[cond] - a) / (b - a))^2
  # Right half: (a+b)/2 < x < b
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 2 * ((x[cond2] - b) / (b - a))^2
  y
}

#' @rdname membership
#' @export
pi_mf = function(x, params) {
  # Compute PI-shaped membership function
  # x: input values, params: c(a, b, c, d) where a < b < c < d
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  y = rep(0, length(x))
  # Plateau
  y[x >= b & x <= C] = 1
  # Rising: left half  a < x <= (a+b)/2
  cond = x > a & x <= (a+b)/2
  y[cond] = 2 * ((x[cond] - a) / (b - a))^2
  # Rising: right half  (a+b)/2 < x < b
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 1 - 2 * ((x[cond2] - b) / (b - a))^2
  # Falling: left half  C < x <= (C+d)/2
  cond3 = x > C & x <= (C+d)/2
  y[cond3] = 1 - 2 * ((x[cond3] - C) / (d - C))^2
  # Falling: right half  (C+d)/2 < x < d
  cond4 = x > (C+d)/2 & x < d
  y[cond4] = 2 * ((x[cond4] - d) / (d - C))^2
  y
}

#' @rdname membership
#' @export
s_mf = function(x, params) {
  # Compute S-shaped membership function
  # x: input values, params: c(a, b) where a < b
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  y = rep(0, length(x))
  y[x >= b] = 1
  # Left half: a < x <= (a+b)/2
  cond = x > a & x <= (a+b)/2
  y[cond] = 2 * ((x[cond] - a) / (b - a))^2
  # Right half: (a+b)/2 < x < b
  cond2 = x > (a+b)/2 & x < b
  y[cond2] = 1 - 2 * ((x[cond2] - b) / (b - a))^2
  y
}

#' @rdname membership
#' @export
plot_mf = function(mf, xlim = c(0, 10), main = NULL) {
  if(!is.function(mf)) stop("mf must be a function.")
  data = data.frame(x = seq(xlim[1], xlim[2], length.out = 200),
                     y = mf(seq(xlim[1], xlim[2], length.out = 200)))
  p = ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1.2) +
    ggplot2::ylim(0, 1.05) + ggplot2::theme_minimal()
  if(!is.null(main)) p = p + ggplot2::labs(title = main)
  p
}
