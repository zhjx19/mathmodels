#' Rank Sum Ratio (RSR) Evaluation
#'
#' @description
#' Performs Rank Sum Ratio (RSR) evaluation on a dataset of positive indicators,
#' computing ranks, weighted RSR values, and a linear regression model to fit RSR
#' against probit-transformed ranks. Supports integer or non-integer ranking methods.
#'
#' @importFrom dplyr mutate across arrange summarise n
#' @importFrom tidyr unnest_longer
#'
#' @param data Data frame with positive indicator data; first column is an ID column
#' for identifying evaluation objects.
#' @param w Numeric vector, weights for indicators (default = equal weights).
#' @param method Character scalar, ranking method: "int" for integer ranks or
#' "non-int" for scaled ranks in [1, n] (default = "int").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{resultTable}: Data frame with RSR values, ranks, cumulative frequencies,
#'   probit values, and fitted RSR values.
#'   \item \code{reg}: Linear model object fitting RSR against probit values.
#'   \item \code{rankTable}: Data frame with ranked indicator values.
#' }
#'
#' @details
#' The \code{rank_sum_ratio} function implements the RSR method for evaluating
#' objects based on positive indicators. It ranks the indicators (using integer or
#' non-integer methods), computes weighted RSR values, adjusts ranks with probit
#' transformation, and fits a linear regression model to relate RSR to probit values.
#' The function assumes the first column of \code{data} is an ID column, and weights
#' (\code{w}) can be provided or set to equal weights by default.
#'
#' @examples
#' # Example data
#' data = data.frame(ID = c("A", "B", "C"), X1 = c(10, 20, 15), X2 = c(5, 10, 8))
#' w = c(0.4, 0.6)
#' rank_sum_ratio(data, w, method = "int")
#'
#' @name rank_sum_ratio
#' @export
rank_sum_ratio = function(data, w = NULL, method = "int") {
  # Compute Rank Sum Ratio (RSR)
  # data: positive indicator data, first column is ID
  # w: weights for indicators
  # method: "int" for integer ranks, "non-int" for scaled ranks
  # Returns: result table, linear regression, rank table
  n = nrow(data)
  m = ncol(data) - 1
  if(is.null(w)) w = rep(1, m)
  # Select ranking method
  switch(method,
         "int" = {
           rankTable = dplyr::mutate(data, dplyr::across(-ID, rank))
         },
         "non-int" = {
           rankTable = dplyr::mutate(data, dplyr::across(-ID, \(x) rescale(x, a = 1, b = n)))
         }
  )
  # Main computation
  rltTable = rankTable |>
    dplyr::mutate(RSR = apply(rankTable[-1], 1, \(x) sum(x * w) / (sum(w) * n)),
                  barR = rank(RSR)) |>
    dplyr::arrange(RSR) |>
    dplyr::summarise(ID = list(ID), f = dplyr::n(), .by = c(RSR, barR)) |>
    dplyr::mutate(sumf = cumsum(f), barRn = barR / n,
                  barRn = ifelse(barRn == 1, 1-1/(4*n), barRn),
                  Probit = 5 + qnorm(barRn))
  reg = lm(RSR ~ Probit, rltTable)
  resultTable = rltTable |>
    dplyr::mutate(RSRfit = predict(reg, rltTable)) |>
    tidyr::unnest_longer(ID) |>
    dplyr::relocate(ID, .before = 1)
  list(resultTable = resultTable, reg = reg, rankTable = rankTable)
}
