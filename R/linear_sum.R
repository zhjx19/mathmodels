#' Linear weighted synthesis
#'
#' @param data Indicator matrix or data frame (rows = samples, columns = indicators)
#' @param w Weight vector (length must match ncol(data))
#'
#' @return Vector of comprehensive scores
#'
#' @export
#'
#' @examples
#' data = data.frame(
#'   GDP = c(0.85, 0.72, 0.91, NA),
#'   Employment = c(0.78, 0.85, 0.67, 0.73),
#'   Environment = c(0.65, 0.72, NA, 0.81)
#' )
#' w = c(0.5, 0.3, 0.2)
#' linear_sum(data, w)

linear_sum = function(data, w) {
  if(!is.data.frame(data) && !is.matrix(data))
    stop("data must be a data frame or matrix.")
  if(!is.numeric(w) || is.matrix(w))
    stop("w must be a numeric vector.")
  if(length(w) != ncol(data))
    stop("Length of w must equal ncol(data).")
  as.matrix(data) %*% w |> as.vector()
}
