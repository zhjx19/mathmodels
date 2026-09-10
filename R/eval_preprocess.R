#' Preprocessing Functions for Data Normalization and Standardization
#'
#' @description
#' A collection of functions to preprocess numeric data, including standardization, L2 norm normalization, Min-Max scaling, centered-type normalization, interval-type normalization, extreme-value-based normalization, initial-value-based normalization, mean-based normalization, and negative-to-positive transformation. These functions transform a numeric vector to a standardized or normalized scale, suitable for various indicator types (positive, negative, centered, interval-based, or extreme-based).
#'
#' @param x Numeric vector to be preprocessed.
#' @param type Character scalar specifying the transformation direction or method:
#' \describe{
#'   \item{"+"}{Positive direction (larger values are better, for `rescale`, `rescale_extreme` and `rescale_initial`).}
#'   \item{"-"}{Negative direction (smaller values are better, for `rescale` `rescale_extreme` and `rescale_initial`).}
#'   \item{"minmax"}{Min-max transformation (for `to_positive`).}
#'   \item{"reciprocal"}{Reciprocal transformation (for `to_positive`).}
#' }
#' @param a Numeric scalar, lower bound of the output range or optimal interval (for `rescale` and `rescale_interval`).
#' @param b Numeric scalar, upper bound of the output range or optimal interval (for `rescale` and `rescale_interval`).
#' @param m Numeric scalar, optimal value for centered-type normalization (for `rescale_middle`).
#' @param center Logical or numeric scalar, passed to `base::scale` for centering (for `standardize`). Default is `TRUE`.
#' @param scale Logical or numeric scalar, passed to `base::scale` for scaling (for `standardize`). Default is `TRUE`.
#'
#' @return A numeric vector of the same length as `x`, transformed as follows:
#' \itemize{
#'   \item `standardize`: Standardized values (mean = 0, sd = 1).
#'   \item `normalize`: L2 norm normalized values (Euclidean norm, unit length).
#'   \item `rescale`: Min-Max scaled values in `[a, b]` (default \code{[0, 1]}).
#'   \item `rescale_middle`: Centered-type normalized values in \code{[0, 1]}, where 1 indicates `x = m`.
#'   \item `rescale_interval`: Interval-type normalized values in \code{[0, 1]}, where 1 indicates `x` in `[a, b]`.
#'   \item `rescale_extreme`: Extreme-based normalized values using `x/max(x)` (positive) or `min(x)/x` (negative).
#'   \item `rescale_initial`: Initial-based normalized values using `x/x[1]` or `x[1]/x`.
#'   \item `rescale_mean`: Mean-based normalized values using `x/mean(x)`.
#'   \item `to_positive`: Converts negative indicators to positive: `minmax` returns `max(x) - x` (a reversal in the original units, not a rescale to \code{[0, 1]}); `reciprocal` returns `1/x` and requires strictly positive data.
#' }
#'
#' @details
#' These functions support various preprocessing needs in data analysis:
#' \itemize{
#'   \item `standardize`: Applies Z-score standardization (mean = 0, sd = 1), ideal for equalizing variances or normally distributed data.
#'   \item `normalize`: Scales the vector to unit length by dividing by its L2 (Euclidean) norm, useful for machine learning or similarity calculations.
#'   \item `rescale`: Performs Min-Max scaling to a specified range (default \code{[0, 1]}), supporting positive or negative indicators.
#'   \item `rescale_middle`: Normalizes centered-type indicators, where values closer to an optimal value `m` are better, mapping to \code{[0, 1]}.
#'   \item `rescale_interval`: Normalizes interval-type indicators, where values within `[a, b]` are optimal, mapping to \code{[0, 1]}.
#'   \item `rescale_extreme`: Normalizes using extreme values: `x/max(x)` for positive indicators or `min(x)/x` for negative indicators, often used in grey relational analysis.
#'   \item `rescale_initial`: Normalizes by dividing by the first value (`x/x[1]` or `x[1]/x`), commonly used in grey relational analysis.
#'   \item `rescale_mean`: Normalizes by dividing by the mean (`x/mean(x)`), commonly used in grey relational analysis.
#'   \item `to_positive`: Converts negative indicators to positive: `minmax` returns `max(x) - x` (a reversal in the original units, not a rescale to `c(0,1)`), and `reciprocal` returns `1/x`, which requires strictly positive data (no zeros or negatives).
#' }
#' Missing values (`NA`) are preserved in the output. For `rescale_initial` and `rescale_mean`, the initial value or mean must be non-zero, respectively.
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
#' rescale(x)                # Scale to \code{[0, 1]}
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
#' rescale_extreme(x)         # x/max(x)
#' rescale_extreme(x, "-")    # min(x)/x
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
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if (isTRUE(scale)) {
    s = stats::sd(x, na.rm = TRUE)
    if (is.na(s)) stop("x has no non-missing values; cannot standardize.")
    if (s == 0) stop("x has zero variance; cannot standardize.")
  }
  as.vector(base::scale(x, center = center, scale = scale))
}

#' @rdname preprocess
#' @export
normalize = function(x) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  nx = norm(x[!is.na(x)], "2")
  if (nx == 0) stop("x has zero L2 norm; cannot normalize.")
  x / nx
}

#' @rdname preprocess
#' @export
rescale = function(x, type = "+", a = 0, b = 1) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!type %in% c("+", "-")) stop('type must be "+" or "-".')
  if (a >= b) stop("a must be < b.")
  rng = range(x, na.rm = TRUE)
  if(rng[2] == rng[1]) stop("x has zero range; cannot rescale.")
  switch(type,
         "+" = (b - a) * (x - rng[1]) / (rng[2] - rng[1]) + a,
         "-" = (b - a) * (rng[2] - x) / (rng[2] - rng[1]) + a)
}

#' @rdname preprocess
#' @export
rescale_middle = function(x, m) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!is.numeric(m) || length(m) != 1) stop("m must be a numeric scalar.")
  # Centered-type rescale, m is the single optimal value
  d = abs(x - m)
  dmax = max(d, na.rm = TRUE)
  if (dmax == 0) return(ifelse(is.na(x), NA_real_, 1))
  1 - d / dmax
}

#' @rdname preprocess
#' @export
rescale_interval = function(x, a, b) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
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
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!type %in% c("+", "-")) stop('type must be "+" or "-".')
  if (any(x <= 0, na.rm = TRUE))
    stop("rescale_extreme() requires strictly positive data (no zeros or negatives).")
  # Normalization using x/max(x) for positive indicators or min(x)/x for negative indicators.
  # Commonly used in grey relational analysis.
  switch(type,
         "+" = x / max(x, na.rm = TRUE),
         "-" = min(x, na.rm = TRUE) / x)
}

#' @rdname preprocess
#' @export
rescale_initial = function(x, type = "+") {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!type %in% c("+", "-")) stop('type must be "+" or "-".')
  # Normalization using x/x[1] for positive indicators, or x[1]/x for negative indicators.
  # Commonly used in grey relational analysis.
  if (is.na(x[1]) || x[1] == 0) stop("x[1] must be non-zero.")
  if (type == "-" & any(x == 0)) stop("x must not contain zeros since reciprocal transformation is used.")
  switch (type,
    "+" = x / x[1],
    "-" = x[1] / x)
}

#' @rdname preprocess
#' @export
rescale_mean = function(x) {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  # Normalization using x/mean(x), commonly used in grey relational analysis.
  mu = mean(x, na.rm = TRUE)
  if (is.na(mu)) stop("x has no non-missing values; mean is NA.")
  if (mu == 0) stop("mean(x) must be non-zero.")
  x / mu
}

#' @rdname preprocess
#' @export
to_positive = function(x, type = "minmax") {
  if(!is.numeric(x)) stop("x must be a numeric vector.")
  if(!type %in% c("minmax", "reciprocal")) stop('type must be "minmax" or "reciprocal".')
  switch(type,
         "minmax" = {
           max(x, na.rm = TRUE) - x
         },
         "reciprocal" = {
           if (any(x <= 0, na.rm = TRUE))
             stop('type = "reciprocal" requires strictly positive data (no zeros or negatives).')
           1 / x
         })
}
