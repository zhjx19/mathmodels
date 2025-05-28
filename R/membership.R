#' @title Membership Functions for Fuzzy Logic
#'
#' @description
#' A collection of functions to compute membership values for various fuzzy sets, including triangular, trapezoidal, Gaussian, generalized bell, two-parameter Gaussian, sigmoid, difference of sigmoids, product of sigmoids, Z-shaped, PI-shaped, and S-shaped membership functions. Includes a function to visualize membership functions using ggplot2. These are designed for evaluation models in mathematical modeling, compatible with \code{fuzzy_eval} in the \code{mathmodels} package.
#'
#' @param x Numeric vector, input values for which to compute membership.
#' @param params Numeric vector, parameters defining the membership function:
#' \itemize{
#'   \item For \code{tri_mf}: \code{c(a, b, c)}, where \code{a <= b <= c} (left base, peak, right base).
#'   \item For \code{trap_mf}: \code{c(a, b, c, d)}, where \code{a <= b <= c <= d} (left base, left top, right top, right base).
#'   \item For \code{gauss_mf}: \code{c(sigma, c)}, where \code{sigma > 0} (spread, center).
#'   \item For \code{gbell_mf}: \code{c(a, b, c)}, where \code{a > 0}, \code{b > 0} (width, shape, center).
#'   \item For \code{gauss2mf}: \code{c(s1, c1, s2, c2)}, where \code{s1 > 0}, \code{s2 > 0} (left spread, left center, right spread, right center).
#'   \item For \code{sigmoid_mf}: \code{c(a, b)}, where \code{a > 0} (slope, inflection point).
#'   \item For \code{dsigmoid_mf}: \code{c(a1, c1, a2, c2)}, where \code{a1 > 0}, \code{a2 > 0} (slopes and inflection points for two sigmoids).
#'   \item For \code{psigmoid_mf}: \code{c(a1, c1, a2, c2)}, where \code{a1 > 0}, \code{a2 > 0} (slopes and inflection points for two sigmoids).
#'   \item For \code{z_mf}: \code{c(a, b)}, where \code{a < b} (left base, right base).
#'   \item For \code{pi_mf}: \code{c(a, b, c, d)}, where \code{a < b < c < d} (left base, left shoulder, right shoulder, right base).
#'   \item For \code{s_mf}: \code{c(a, b)}, where \code{a < b} (left base, right base).
#' }
#' @param mf Function, a membership function with fixed parameters (e.g., \code{function(x) tri_mf(x, c(2, 5, 8))}).
#' @param xlim Numeric vector of length 2, x-axis limits for plotting (default \code{c(0, 10)}).
#' @param main Character, plot title (default \code{NULL}, no title).
#'
#' @return
#' \itemize{
#'   \item For membership functions (\code{tri_mf}, \code{trap_mf}, \code{gauss_mf}, \code{gbell_mf}, \code{gauss2mf}, \code{sigmoid_mf}, \code{dsigmoid_mf}, \code{psigmoid_mf}, \code{z_mf}, \code{pi_mf}, \code{s_mf}): A numeric vector of membership values in [0, 1], same length as \code{x}.
#'   \item For \code{plot_mf}: A ggplot2 object, plotting the membership function.
#' }
#'
#' @details
#' These functions support evaluation models in mathematical modeling:
#' \itemize{
#'   \item \code{tri_mf}: Triangular membership, linear rise from \code{a} to \code{b} (peak) and fall to \code{c}.
#'   \item \code{trap_mf}: Trapezoidal membership, linear rise from \code{a} to \code{b}, plateau from \code{b} to \code{c}, fall to \code{d}.
#'   \item \code{gauss_mf}: Gaussian membership, bell-shaped curve centered at \code{c} with spread \code{sigma}.
#'   \item \code{gbell_mf}: Generalized bell membership, bell-shaped curve with width \code{a}, shape \code{b}, and center \code{c}.
#'   \item \code{gauss2mf}: Two-parameter Gaussian membership, combining two Gaussians with spreads \code{s1}, \code{s2} and centers \code{c1}, \code{c2}.
#'   \item \code{sigmoid_mf}: Sigmoid membership, S-shaped curve with slope \code{a} and inflection point \code{b}.
#'   \item \code{dsigmoid_mf}: Difference of two sigmoids, combining slopes \code{a1}, \code{a2} and inflection points \code{c1}, \code{c2}.
#'   \item \code{psigmoid_mf}: Product of two sigmoids, combining slopes \code{a1}, \code{a2} and inflection points \code{c1}, \code{c2}.
#'   \item \code{z_mf}: Z-shaped membership, decreasing from 1 at \code{a} to 0 at \code{b}.
#'   \item \code{pi_mf}: PI-shaped membership, rising from \code{a} to \code{b}, plateau from \code{b} to \code{c}, falling to \code{d}.
#'   \item \code{s_mf}: S-shaped membership, increasing from 0 at \code{a} to 1 at \code{b}.
#'   \item \code{plot_mf}: Plots a membership function over \code{xlim} using ggplot2, suitable for tidyverse workflows.
#' }
#' Membership values can be used to construct fuzzy evaluation matrices for \code{fuzzy_eval}.
#' Implemented in base R, except \code{plot_mf}, which requires ggplot2.
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
  sigma = params[1]
  C = params[2]
  exp(-(x - C)^2 / (2 * sigma^2))
}

#' @rdname membership
#' @export
gbell_mf = function(x, params) {
  # Compute generalized bell membership function
  # x: input values, params: c(a, b, c) where a > 0, b > 0, c > 0
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
  a1 = params[1]
  c1 = params[2]
  sigmoid(x, a1, c1)
}


#' @rdname membership
#' @export
dsigmoid_mf = function(x, params) {
  # Compute difference of sigmoid membership function
  # x: input values, params: c(a1, c1, a2, c2) where a1 > 0, a2 > 0
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
  a = params[1]
  b = params[2]
  y = rep(1, length(x))
  idx1 = a<=x & x<=(a+b)/2
  idx2 = (a+b)/2<=x & x<=b
  y[idx1] = 1- 2 * ((x[idx1]-a)/(b-a))^2
  y[idx2] = 2 * ((x[idx2]-b)/(b-a))^2
  y[x >= b] = 0
  y
}


#' @rdname membership
#' @export
pi_mf = function(x, params) {
  # Compute PI-shaped membership function
  # x: input values, params: c(a, b, C, d) where a < b < C < d
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  y = rep(0, length(x))
  idx1 = a<=x & x<=(a+b)/2
  idx2 = (a+b)/2<=x & x<=b
  idx3 = C<=x & x<=(C+d)/2
  idx4 = (C+d)/2<=x & x<=d
  y[idx1] = 2 * ((x[idx1]-a)/(b-a))^2
  y[idx2] = 1 - 2 * ((x[idx2]-b)/(b-a))^2
  y[b<=x & x<=C] = 1
  y[idx3] = 1 - 2 * ((x[idx3]-C)/(d-C))^2
  y[idx4] = 2 * ((x[idx4]-d)/(d-C))^2
  y
}

#' @rdname membership
#' @export
s_mf = function(x, params) {
  # Compute S-shaped membership function
  # x: input values, params: c(a, b) where a < b
  a = params[1]
  b = params[2]
  y = rep(0, length(x))
  idx1 = a<=x & x<=(a+b)/2
  idx2 = (a+b)/2<=x & x<=b
  y[idx1] = 2 * ((x[idx1]-a)/(b-a))^2
  y[idx2] = 1 - 2 * ((x[idx2]-b)/(b-a))^2
  y[x>=b] = 1
  y
}


#' @rdname membership
#' @export
plot_mf = function(mf, xlim = c(0, 10), main = NULL) {
  # Plot a membership function
  # mf: membership function, xlim: x-axis limits, main: plot title
  ggplot2::ggplot() +
    ggplot2::geom_function(fun = mf, xlim = xlim,
                           color = "red", linewidth = 1.2) +
    ggplot2::scale_x_continuous(breaks = pretty(xlim)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = main) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
