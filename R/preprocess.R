#' Preprocessing Functions for Data Normalization and Standardization
#'
#' @description
#' A collection of functions for preprocessing numeric data, including standardization,
#' L2 norm normalization, Min-Max normalization, centered-type normalization,
#' interval-type normalization, and negative-to-positive transformation. Each function
#' transforms a numeric vector to a standardized or normalized scale, handling different
#' types of indicators (positive, negative, centered, or interval-based).
#'
#' @param X Numeric vector to be preprocessed.
#' @param type Character scalar indicating the transformation direction or type:
#' \describe{
#'   \item{"+"}{Positive direction (larger values are better, for \code{rescale}).}
#'   \item{"-"}{Negative direction (smaller values are better, for \code{rescale}).}
#'   \item{"reciprocal"}{Reciprocal transformation (for \code{neg_to_pos}).}
#'   \item{"minmax"}{Min-max transformation (for \code{neg_to_pos}).}
#' }
#' @param switch Character scalar indicating the specific transformation for \code{neg_to_pos}:
#' \describe{
#'   \item{"reciprocal"}{Applies reciprocal transformation (1/x).}
#'   \item{"minmax"}{Applies min-max transformation (max(x) - x).}
#' }
#' @param a Numeric scalar, lower bound of the output range or interval (for \code{rescale}
#' and \code{rescale_interval}).
#' @param b Numeric scalar, upper bound of the output range or interval (for \code{rescale}
#' and \code{rescale_interval}).
#' @param m Numeric scalar, the optimal value for centered-type normalization (for \code{rescale_middle}).
#' @param center Logical or numeric scalar, passed to \code{base::scale} for centering
#' (for \code{standardize}). Default is \code{TRUE}.
#' @param scale Logical or numeric scalar, passed to \code{base::scale} for scaling
#' (for \code{standardize}). Default is \code{TRUE}.
#'
#' @return A numeric vector of the same length as \code{x}, transformed according to
#' the specified method:
#' \itemize{
#'   \item \code{standardize}: Standardized values (mean = 0, sd = 1).
#'   \item \code{normalize}: Normalized values using L2 norm (Euclidean norm).
#'   \item \code{rescale}: Normalized values in \code{[a, b]} (default [0, 1]).
#'   \item \code{rescale_middle}: Normalized values in [0, 1], where 1 indicates \code{x = m}.
#'   \item \code{rescale_interval}: Normalized values in [0, 1], where 1 indicates \code{x} in \code{[a, b]}.
#'   \item \code{to_positive}: Transformed values where negative indicators are converted to positive
#'   using either reciprocal or min-max transformation.
#' }
#'
#' @details
#' These functions are tailored for different indicator types in data analysis:
#' \itemize{
#'   \item \code{standardize}: Applies Z-score standardization, transforming data
#'   to have mean 0 and standard deviation 1. Suitable for normally distributed data or
#'   when equalizing variances.
#'   \item \code{normalize}: Normalizes data by dividing by the L2 (Euclidean) norm, scaling
#'   the vector to unit length. Useful for machine learning or similarity computations.
#'   \item \code{rescale}: Performs Min-Max normalization, scaling data to a specified range
#'   (default [0, 1]). Supports positive or negative indicators.
#'   \item \code{rescale_middle}: Normalizes centered-type indicators, where values closer to
#'   an optimal value \code{m} are better. Output is in [0, 1].
#'   \item \code{rescale_interval}: Normalizes interval-type indicators, where values in the
#'   optimal interval \code{[a, b]} are best. Output is in [0, 1].
#'   \item \code{to_positive}: Converts negative indicators to positive using either reciprocal
#'   transformation (1/x) or min-max transformation (max(x) - x). The \code{type} and \code{switch}
#'   parameters must match (e.g., both "reciprocal" or both "minmax").
#' }
#'
#' @examples
#' # Standardization
#' x = c(4, 1, NA, 5, 8)
#' standardize(x)
#'
#' # L2 norm normalization
#' normalize(x)
#'
#' # Min-Max normalization (positive direction)
#' rescale(x)                # Scale to [0, 1]
#' rescale(x, type = "-", a = 0.002, b = 0.996)  # Reverse scaling
#'
#' # Negative-to-positive transformation
#' to_positive(x)                       # Min-max transformation
#' to_positive(x, type = "reciprocal")  # Reciprocal transformation
#'
#' # Centered-type normalization
#' PH = 6:9
#' rescale_middle(PH, 7)
#'
#' # Interval-type normalization
#' Temp = c(35.2, 35.8, 36.6, 37.1, 37.8, 38.4)
#' rescale_interval(Temp, 36, 37)
#'
#' @name preprocess
NULL

#' @rdname preprocess
#' @export
standardize = function(x, center = TRUE, scale = TRUE) {
  as.vector(base::scale(x, center = center, scale = scale))
}

#' @rdname preprocess
#' @export
normalize = function(x) {
  x / norm(x[!is.na(x)], "2")
}

#' @rdname preprocess
#' @export
rescale = function(x, type = "+", a = 0, b = 1) {
  rng = range(x, na.rm = TRUE)
  switch(type,
         "+" = (b - a) * (x - rng[1]) / (rng[2] - rng[1]) + a,
         "-" = (b - a) * (rng[2] - x) / (rng[2] - rng[1]) + a)
}

#' @rdname preprocess
#' @export
rescale_middle = function(x, m) {
  # Centered-type rescale, m is the single optimal value
  d = abs(x - m)
  1 - d / max(d, na.rm = TRUE)
}

#' @rdname preprocess
#' @export
rescale_interval = function(x, a, b) {
  # Interval-type rescale, [a,b] is the optimal interval
  M = max(a - min(x, na.rm=TRUE), max(x, na.rm=TRUE) - b)
  y = rep(1, length(x))
  y[x<a & !is.na(x)] = 1 - (a-x[x<a & !is.na(x)]) / M
  y[x>b & !is.na(x)] = 1 - (x[x>b & !is.na(x)]-b) / M
  y[is.na(x)] = NA
  y
}

#' @rdname preprocess
#' @export
to_positive = function(x, type = "minmax") {
  switch(type,
         "minmax" = {
           max(x, na.rm = TRUE) - x
         },
         "reciprocal" = {
           1 / x
         })
}
