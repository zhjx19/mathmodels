#' Preprocessing Functions for Data Normalization and Standardization
#'
#' @description
#' A collection of functions to preprocess numeric data, including standardization, L2 norm normalization, Min-Max scaling, centered-type normalization, interval-type normalization, extreme-value-based normalization, initial-value-based normalization, mean-based normalization, and negative-to-positive transformation. These functions transform a numeric vector to a standardized or normalized scale, suitable for various indicator types (positive, negative, centered, interval-based, or extreme-based).
#'
#' @param x Numeric vector to be preprocessed.
#' @param type Character scalar specifying the transformation direction or method:
#' \describe{
#'   \item{"+"}{Positive direction (larger values are better, for \code{rescale}, \code{rescale_extreme} and \code{rescale_initial}).}
#'   \item{"-"}{Negative direction (smaller values are better, for \code{rescale} \code{rescale_extreme} and \code{rescale_initial}).}
#'   \item{"minmax"}{Min-max transformation (for \code{to_positive}).}
#'   \item{"reciprocal"}{Reciprocal transformation (for \code{to_positive}).}
#' }
#' @param a Numeric scalar, lower bound of the output range or optimal interval (for \code{rescale} and \code{rescale_interval}).
#' @param b Numeric scalar, upper bound of the output range or optimal interval (for \code{rescale} and \code{rescale_interval}).
#' @param m Numeric scalar, optimal value for centered-type normalization (for \code{rescale_middle}).
#' @param center Logical or numeric scalar, passed to \code{base::scale} for centering (for \code{standardize}). Default is \code{TRUE}.
#' @param scale Logical or numeric scalar, passed to \code{base::scale} for scaling (for \code{standardize}). Default is \code{TRUE}.
#'
#' @return A numeric vector of the same length as \code{x}, transformed as follows:
#' \itemize{
#'   \item \code{standardize}: Standardized values (mean = 0, sd = 1).
#'   \item \code{normalize}: L2 norm normalized values (Euclidean norm, unit length).
#'   \item \code{rescale}: Min-Max scaled values in \code{[a, b]} (default [0, 1]).
#'   \item \code{rescale_middle}: Centered-type normalized values in [0, 1], where 1 indicates \code{x = m}.
#'   \item \code{rescale_interval}: Interval-type normalized values in [0, 1], where 1 indicates \code{x} in \code{[a, b]}.
#'   \item \code{rescale_extreme}: Extreme-based normalized values using \code{min(x)/x} (positive) or \code{x/max(x)} (negative).
#'   \item \code{rescale_initial}: Initial-based normalized values using \code{x/x[1]} or \code{x[1]/x}.
#'   \item \code{rescale_mean}: Mean-based normalized values using \code{x/mean(x)}.
#'   \item \code{to_positive}: Transformed values converting negative indicators to positive using min-max or reciprocal transformation.
#' }
#'
#' @details
#' These functions support various preprocessing needs in data analysis:
#' \itemize{
#'   \item \code{standardize}: Applies Z-score standardization (mean = 0, sd = 1), ideal for equalizing variances or normally distributed data.
#'   \item \code{normalize}: Scales the vector to unit length by dividing by its L2 (Euclidean) norm, useful for machine learning or similarity calculations.
#'   \item \code{rescale}: Performs Min-Max scaling to a specified range (default [0, 1]), supporting positive or negative indicators.
#'   \item \code{rescale_middle}: Normalizes centered-type indicators, where values closer to an optimal value \code{m} are better, mapping to [0, 1].
#'   \item \code{rescale_interval}: Normalizes interval-type indicators, where values within \code{[a, b]} are optimal, mapping to [0, 1].
#'   \item \code{rescale_extreme}: Normalizes using extreme values: \code{min(x)/x} for positive indicators or \code{x/max(x)} for negative indicators, often used in grey relational analysis.
#'   \item \code{rescale_initial}: Normalizes by dividing by the first value (\code{x/x[1]} or \code{x[1]/x}), commonly used in grey relational analysis.
#'   \item \code{rescale_mean}: Normalizes by dividing by the mean (\code{x/mean(x)}), commonly used in grey relational analysis.
#'   \item \code{to_positive}: Converts negative indicators to positive using either min-max (\code{max(x) - x}) or reciprocal (\code{1/x}) transformation.
#' }
#' Missing values (\code{NA}) are preserved in the output. For \code{rescale_initial} and \code{rescale_mean}, the initial value or mean must be non-zero, respectively.
#'
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
#' # Extreme-based normalization
#' rescale_extreme(x)         # min(x)/x
#' rescale_extreme(x, "-")    # x/max(x)
#'
#' # Initial-based normalization
#' rescale_initial(x)
#'
#' # Mean-based normalization
#' rescale_mean(x)
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
rescale_extreme = function(x, type = "+") {
  # Normalization using x/max(x) for positive indicators or min(x)/x for negative indicators.
  # Commonly used in grey relational analysis.
  switch(type,
         "+" = x / max(x, na.rm = TRUE),
         "-" = min(x, na.rm = TRUE) / x)
}

#' @rdname preprocess
#' @export
rescale_initial = function(x, type = "+") {
  # Normalization using x/x[1] for positive indicators, or x[1]/x for negative indicators.
  # Commonly used in grey relational analysis.
  if (x[1] == 0 | is.na(x[1])) stop("x[1] must be non-zero.")
  if (type == "-" & any(x == 0)) stop("x must not contain zeros since reciprocal transformation is used.")
  switch (type,
    "+" = x / x[1],
    "-" = x[1] / x)
}

#' @rdname preprocess
#' @export
rescale_mean = function(x) {
  # Normalization using x/mean(x), commonly used in grey relational analysis.
  mu = mean(x, na.rm = TRUE)
  if (mu == 0) stop("mean(x) must be non-zero.")
  x / mu
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
