# Evaluation Algorithms — Fuzzy Evaluation
#
# Includes: membership functions, fuzzy_eval, compute_mf, and defuzzify
#


# ============================================================================
# Membership Functions
# ============================================================================

#' @title Membership Functions for Fuzzy Logic
#'
#' @description
#' A collection of functions to compute membership values for various fuzzy sets, including triangular, trapezoidal, Gaussian, generalized bell, two-parameter Gaussian, sigmoid, difference of sigmoids, product of sigmoids, Z-shaped, PI-shaped, and S-shaped membership functions. Includes a function to visualize membership functions using ggplot2. These are designed for evaluation models in mathematical modeling, compatible with `fuzzy_eval` in the `mathmodels` package.
#'
#' @param x Numeric vector, input values for which to compute membership.
#' @param params Numeric vector, parameters defining the membership function:
#' \itemize{
#'   \item For `mf_tri`: `c(a, b, c)`, where `a < b < c` (left base, peak, right base).
#'   \item For `mf_trap`: `c(a, b, c, d)`, where `a < b <= c < d` (left base, left top,
#'     right top, right base); one-sided shoulders are supported with
#'     `a = b = -Inf` (no rise) or `c = d = Inf` (no fall).
#'   \item For `mf_gauss`: `c(sigma, c)`, where `sigma > 0` (spread, center).
#'   \item For `mf_gbell`: `c(a, b, c)`, where `a > 0`, `b > 0` (width, shape, center).
#'   \item For `mf_gauss2`: `c(s1, c1, s2, c2)`, where `s1 > 0`, `s2 > 0` (left spread, left center, right spread, right center).
#'   \item For `mf_sigmoid`: `c(a, b)`, where `a > 0` (slope, inflection point).
#'   \item For `mf_dsigmoid`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0` (slopes and inflection points for two sigmoids).
#'   \item For `mf_psigmoid`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0` (slopes and inflection points for two sigmoids).
#'   \item For `mf_z`: `c(a, b)`, where `a < b` (left base, right base).
#'   \item For `mf_pi`: `c(a, b, c, d)`, where `a < b <= c < d` (left base, left shoulder, right shoulder, right base).
#'   \item For `mf_s`: `c(a, b)`, where `a < b` (left base, right base).
#' }
#' @param mf Function, a membership function with fixed parameters (e.g., `function(x) mf_tri(x, c(2, 5, 8))`).
#' @param xlim Numeric vector of length 2, x-axis limits for plotting (default `c(0, 10)`).
#' @param main Character, plot title (default `NULL`, no title).
#'
#' @return
#' \itemize{
#'   \item For membership functions (`mf_tri`, `mf_trap`, `mf_gauss`, `mf_gbell`, `mf_gauss2`, `mf_sigmoid`, `mf_dsigmoid`, `mf_psigmoid`, `mf_z`, `mf_pi`, `mf_s`): A numeric vector of membership values in \code{[0, 1]}, same length as `x`.
#'   \item For `plot_mf`: A ggplot2 object, plotting the membership function.
#' }
#'
#' @details
#' These functions support evaluation models in mathematical modeling:
#' \itemize{
#'   \item `mf_tri`: Triangular membership, linear rise from `a` to `b` (peak) and fall to `c`.
#'   \item `mf_trap`: Trapezoidal membership, linear rise from `a` to `b`, plateau from `b` to `c`, fall to `d`.
#'   \item `mf_gauss`: Gaussian membership, bell-shaped curve centered at `c` with spread `sigma`.
#'   \item `mf_gbell`: Generalized bell membership, bell-shaped curve with width `a`, shape `b`, and center `c`.
#'   \item `mf_gauss2`: Two-parameter Gaussian membership, combining two Gaussians with spreads `s1`, `s2` and centers `c1`, `c2`.
#'   \item `mf_sigmoid`: Sigmoid membership, S-shaped curve with slope `a` and inflection point `b`.
#'   \item `mf_dsigmoid`: Difference of two sigmoids, combining slopes `a1`, `a2` and inflection points `c1`, `c2`.
#'   \item `mf_psigmoid`: Product of two sigmoids, combining slopes `a1`, `a2` and inflection points `c1`, `c2`.
#'   \item `mf_z`: Z-shaped membership, decreasing from 1 at `a` to 0 at `b`.
#'   \item `mf_pi`: PI-shaped membership, rising from `a` to `b`, plateau from `b` to `c`, falling to `d`.
#'   \item `mf_s`: S-shaped membership, increasing from 0 at `a` to 1 at `b`.
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
#' mf_tri(x, params = c(3, 6, 8))
#'
#' # Trapezoidal membership
#' mf_trap(x, params = c(1, 5, 7, 8))
#'
#' # Gaussian membership
#' mf_gauss(x, params = c(2, 5))
#'
#' # Generalized bell membership
#' mf_gbell(x, params = c(2, 4, 6))
#'
#' # Two-parameter Gaussian membership
#' mf_gauss2(x, params = c(1, 3, 3, 4))
#'
#' # Sigmoid membership
#' mf_sigmoid(x, params = c(2, 4))
#'
#' # Difference of sigmoids membership
#' dsigmoid_y = mf_dsigmoid(x, params = c(5, 2, 5, 7))
#'
#' # Product of sigmoids membership
#' mf_psigmoid(x, params = c(2, 3, -5, 8))
#'
#' # Z-shaped membership
#' mf_z(x, params = c(3, 7))
#'
#' # PI-shaped membership
#' mf_pi(x, params = c(1, 4, 5, 10))
#'
#' # S-shaped membership
#' mf_s(x, params = c(1, 8))
#'
#' \dontrun{
#' # Visualize membership functions
#' plot_mf(\(x) mf_tri(x, c(3, 6, 8)), main = "Triangular MF")
#' plot_mf(\(x) mf_trap(x, c(1, 5, 7, 8)), main = "Trapezoidal MF")
#' plot_mf(\(x) mf_gauss(x, c(2, 5)), main = "Gaussian MF")
#' plot_mf(\(x) mf_gbell(x, c(2, 4, 6)), main = "Generalized Bell MF")
#' plot_mf(\(x) mf_gauss2(x, c(1, 3, 3, 4)), main = "Two-Parameter Gaussian MF")
#' plot_mf(\(x) mf_sigmoid(x, c(2, 4)), main = "Sigmoid MF")
#' plot_mf(\(x) mf_dsigmoid(x, c(5, 2, 5, 7)), main = "Difference of Sigmoids MF")
#' plot_mf(\(x) mf_psigmoid(x, c(2, 3, -5, 8)), main = "Product of Sigmoids MF")
#' plot_mf(\(x) mf_z(x, c(3, 7)), main = "Z-Shaped MF")
#' plot_mf(\(x) mf_pi(x, c(1, 4, 5, 10)), main = "PI-Shaped MF")
#' plot_mf(\(x) mf_s(x, c(1, 8)), main = "S-Shaped MF")
#' }
#'
#' @name membership
NULL

#' @rdname membership
#' @export
mf_tri = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 3)
    stop("params must be a numeric vector of length 3.")
  a = params[1]
  b = params[2]
  C = params[3]
  if(!(a < b && b < C))
    stop("params must satisfy a < b < c.", call. = FALSE)
  y = rep(0, length(x))
  y[x >= a & x <= b] = (x[x >= a & x <= b] - a) / (b - a)
  y[x > b & x <= C] = (C - x[x > b & x <= C]) / (C - b)
  y
}

#' @rdname membership
#' @export
mf_trap = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  if(!(a < b || (a == b && is.infinite(a) && is.infinite(b))) ||
     !(b <= C) ||
     !(C < d || (C == d && is.infinite(C) && is.infinite(d))))
    stop("params must satisfy a < b <= c < d.", call. = FALSE)
  y = rep(0, length(x))
  y[x >= a & x <= b] = (x[x >= a & x <= b] - a) / (b - a)
  y[x > b & x < C] = 1
  y[x >= C & x <= d] = (d - x[x >= C & x <= d]) / (d - C)
  y
}

#' @rdname membership
#' @export
mf_gauss = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  sigma = params[1]
  C = params[2]
  if(!(sigma > 0))
    stop("sigma must be positive.", call. = FALSE)
  exp(-(x - C)^2 / (2 * sigma^2))
}

#' @rdname membership
#' @export
mf_gbell = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 3)
    stop("params must be a numeric vector of length 3.")
  a = params[1]
  b = params[2]
  C = params[3]
  if(!(a > 0))
    stop("a must be positive.", call. = FALSE)
  1 / (1 + abs((x - C) / a)^(2*b))
}

#' @rdname membership
#' @export
mf_gauss2 = function(x, params) {
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
mf_sigmoid = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a1 = params[1]
  c1 = params[2]
  sigmoid(x, a1, c1)
}

#' @rdname membership
#' @export
mf_dsigmoid = function(x, params) {
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
mf_psigmoid = function(x, params) {
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
mf_z = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  if(!(a < b))
    stop("params must satisfy a < b.", call. = FALSE)
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
mf_pi = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 4)
    stop("params must be a numeric vector of length 4.")
  a = params[1]
  b = params[2]
  C = params[3]
  d = params[4]
  if(!(a < b && b <= C && C < d))
    stop("params must satisfy a < b <= c < d.", call. = FALSE)
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
mf_s = function(x, params) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(params) || length(params) != 2)
    stop("params must be a numeric vector of length 2.")
  a = params[1]
  b = params[2]
  if(!(a < b))
    stop("params must satisfy a < b.", call. = FALSE)
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
# compute_mf — indicator value to membership vector
# ============================================================================

# ---------- internal builder: triangular / trapezoidal (piecewise linear) ----------

build_mf_tri = function(knots, ...) {
  n = length(knots)
  mfs = vector("list", n)

  mfs[[1]] = function(x) mf_trap(x, c(-Inf, -Inf, knots[1], knots[2]))

  mfs[[n]] = function(x) mf_trap(x, c(knots[n - 1], knots[n], Inf, Inf))

  if (n > 2) {
    i = seq_len(n - 2)
    mfs[2:(n-1)] = lapply(i, function(i) {
      t = knots[i + 0:2]
      function(x) mf_tri(x, t)
    })
  }

  mfs
}

# ---------- internal builder: Gaussian ----------

build_mf_gauss = function(knots, sigma = NULL, ...) {
  n = length(knots)
  if (is.null(sigma)) sigma = mean(diff(knots)) / 3
  if (!is.numeric(sigma) || length(sigma) != 1 || sigma <= 0)
    stop("sigma must be a positive numeric scalar.")

  mfs = vector("list", n)

  c1 = knots[1]
  mfs[[1]] = function(x) {
    y = rep(1, length(x))
    y[x > c1] = exp(-(x[x > c1] - c1)^2 / (2 * sigma^2))
    y
  }

  if (n > 2) {
    i = seq_len(n - 2)
    mfs[2:(n-1)] = lapply(i, function(i) {
      function(x) mf_gauss(x, c(sigma, knots[i + 1]))
    })
  }

  cn = knots[n]
  mfs[[n]] = function(x) {
    y = exp(-(x - cn)^2 / (2 * sigma^2))
    y[x >= cn] = 1
    y
  }

  mfs
}

# ---------- internal builder: Sigmoid ----------

build_mf_sigmoid = function(knots, slope = NULL, ...) {
  n = length(knots)
  if (is.null(slope)) slope = 5 / mean(diff(knots))
  if (!is.numeric(slope) || length(slope) != 1 || slope <= 0)
    stop("slope must be a positive numeric scalar.")

  mfs = vector("list", n)

  mfs[[1]] = function(x) mf_sigmoid(x, c(-slope, knots[1]))

  if (n > 2) {
    i = seq_len(n - 2)
    mfs[2:(n-1)] = lapply(i, function(i) {
      cen = knots[i + 1]
      delta = min(cen - knots[i], knots[i + 2] - cen) / 3
      function(x) mf_dsigmoid(x, c(slope, cen - delta, slope, cen + delta))
    })
  }

  mfs[[n]] = function(x) mf_sigmoid(x, c(slope, knots[n]))

  mfs
}

# ---------- public API ----------

#' Build membership functions from knots
#'
#' @description
#' `compute_mf_funs` constructs a list of membership functions (one per evaluation
#' level) from knot values. `compute_mf` evaluates a single indicator value
#' against those functions, returning a membership vector.
#'
#' @inheritParams membership
#' @param knots A numeric vector of length \eqn{n \ge 2} defining the
#'   characteristic values (peaks / centers) of each evaluation level.
#' @param .builder A character string or function that specifies how membership
#'   functions are built from `knots`:
#'   \describe{
#'     \item{`"tri"`}{(default) Piecewise linear (triangular + trapezoidal).
#'       First level decays linearly from 1 to 0 across \code{[th[1], th[2]]},
#'       last level rises from 0 to 1 across \code{[th[n-1], th[n]]}, middle
#'       levels are isosceles triangles peaking at each knot.}
#'     \item{`"gauss"`}{Gaussian (normal) membership. First level is a
#'       right-half Gaussian decaying from 1 at \code{th[1]}, last level is a
#'       left-half Gaussian rising to 1 at \code{th[n]}, middle levels are full
#'       Gaussians centered at each knot.  Pass \code{sigma} to control
#'       spread (default = \code{mean(diff(knots)) / 3}).}
#'     \item{`"sigmoid"`}{Sigmoid-based membership. First level uses a
#'       decreasing sigmoid, last level an increasing sigmoid, middle levels
#'       use difference-of-sigmoids (bell-shaped). Pass \code{slope} to control
#'       steepness (default = \code{5 / mean(diff(knots))}).}
#'     \item{User-supplied function}{A function with signature
#'       \code{function(knots, ...)} that returns a list of \eqn{n}
#'       functions, each accepting a numeric vector \code{x} and returning
#'       membership values in \eqn{[0, 1]}.}
#'   }
#' @param ... Additional arguments passed to the builder (e.g. \code{sigma}
#'   for \code{"gauss"}, \code{slope} for \code{"sigmoid"}, or forwarded to a
#'   custom builder function).
#'
#' @return
#' \itemize{
#'   \item `compute_mf_funs`: A list of \eqn{n} functions, one per level.
#'   \item `compute_mf`: A numeric vector of length \eqn{n} with membership
#'     degrees for each level.
#' }
#'
#' @examples
#' # Triangular membership (default)
#' th = c(0.05, 0.15, 0.25, 0.5)
#' compute_mf(0.07, th)
#'
#' # Gaussian membership with custom sigma
#' compute_mf(0.07, th, .builder = "gauss", sigma = 0.05)
#'
#' # Sigmoid membership with custom slope
#' compute_mf(0.07, th, .builder = "sigmoid", slope = 20)
#'
#' # Custom builder: exponential decay
#' exp_builder = function(knots, rate = 1) {
#'   n = length(knots)
#'   lapply(seq_len(n), function(i) {
#'     force(i)
#'     function(x) exp(-rate * abs(x - knots[i]))
#'   })
#' }
#' compute_mf(0.07, th, .builder = exp_builder, rate = 10)
#'
#' \dontrun{
#' # Visualise all levels for a given builder
#' mfs = compute_mf_funs(th, .builder = "gauss", sigma = 0.05)
#' plots = lapply(mfs, \(f) plot_mf(f, xlim = c(0, 0.6)))
#' gridExtra::grid.arrange(grobs = plots, nrow = 2)
#' }
#'
#' @name compute_mf
NULL

#' @rdname compute_mf
#' @export
compute_mf_funs = function(knots, .builder = "tri", ...) {
  if (!is.numeric(knots) || is.matrix(knots))
    stop("knots must be a numeric vector.")
  n = length(knots)
  if (n < 2) stop("knots must contain at least two values")

  if (is.function(.builder)) {
    return(.builder(knots, ...))
  }

  if (!is.character(.builder) || length(.builder) != 1)
    stop(".builder must be a character string or a function.")

  .builder = match.arg(.builder, c("tri", "gauss", "sigmoid"))

  mfs = switch(.builder,
    tri     = build_mf_tri(knots, ...),
    gauss   = build_mf_gauss(knots, ...),
    sigmoid = build_mf_sigmoid(knots, ...)
  )

  mfs
}

#' @rdname compute_mf
#' @export
compute_mf = function(x, knots, .builder = "tri", ...) {
  if (!is.numeric(x) || length(x) != 1)
    stop("x must be a numeric scalar.")
  if (!is.numeric(knots) || is.matrix(knots))
    stop("knots must be a numeric vector.")
  n = length(knots)
  if (n < 2) stop("knots must contain at least two values")

  mfs = compute_mf_funs(knots, .builder = .builder, ...)
  sapply(mfs, function(f) f(x))
}


# ============================================================================
# fuzzy_eval — fuzzy comprehensive evaluation
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
#' as evaluation grades. Values must be in \code{[0, 1]} with no NA.
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
#' Note: type 5 is a weighted-sum variant (a normalized min-sum); it is not one
#' of the classic Zadeh (min-max) operator pairs.
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
  if(!is.matrix(R) || !is.numeric(R))
    stop("R must be a numeric matrix.", call. = FALSE)
  if(anyNA(R) || any(R < 0) || any(R > 1))
    stop("R must contain only values between 0 and 1, with no NA.", call. = FALSE)
  if(length(w) != ncol(R))
    stop("Length of w must equal ncol(R).")
  if(anyNA(w) || any(w < 0))
    stop("w must be non-negative and contain no NA.", call. = FALSE)
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
# defuzzify — defuzzification
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
#' @details
#' For `max_membership`, when several levels tie for the maximum membership,
#' the score of the first (lowest-index) maximum is returned.
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
         max_membership = scores[which.max(mu)],
         centroid = sum(mu * scores) / sum(mu),
         stop('method must be one of "weighted_average", "max_membership", "centroid"',
              call. = FALSE)
  )
}
