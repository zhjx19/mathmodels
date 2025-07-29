#' Regional Economics Functions
#'
#' @description
#' A collection of functions for calculating regional economics indices, including
#' Location Quotient (LQ), Herfindahl-Hirschman Index (HHI), and Ellison-Glaeser Index (EG).
#' These functions are designed for analyzing regional and industrial data to assess
#' spatial concentration and market structure.
#'
#' @importFrom dplyr mutate across arrange summarise reframe
#' @importFrom tidyr pivot_wider
#' @importFrom tibble enframe
#' @importFrom purrr map2_dbl
#'
#' @param data A data frame containing the necessary data for calculations. For \code{LQ},
#' the first row is assumed to be the national total, with the first two columns specifying
#' the region and total, optionally including a grouping column. For \code{EG}, it contains
#' region, industry, and indicator columns (e.g., employment or output).
#' @param region The name of the column in \code{data} specifying the region (used in \code{LQ} and \code{EG}).
#' @param total The name of the column in \code{data} specifying the total (used in \code{LQ}).
#' @param .cols The columns in \code{data} for which to calculate the location quotient (used in \code{LQ}).
#' @param .by Optional grouping column for \code{LQ}, defaults to \code{NULL} (no grouping).
#' @param x A numeric vector for calculating the HHI (used in \code{HHI}).
#' @param scaled Logical; if \code{TRUE}, the HHI is scaled to account for the number of firms. Defaults to \code{FALSE}.
#' @param industry The name of the column in \code{data} specifying the industry (used in \code{EG}).
#' @param y The name of the column in \code{data} specifying the indicator (e.g., employment or output, used in \code{EG}).
#'
#' @details
#' \describe{
#'   \item{LQ}{Calculates the Location Quotient for multiple columns with optional grouping.
#'     The LQ measures the relative concentration of an industry in a region compared to a national benchmark.
#'     The function assumes the first row of \code{data} contains national totals, with the first two columns
#'     specifying the region and total, and the remaining columns used for LQ calculation.}
#'   \item{HHI}{Calculates the Herfindahl-Hirschman Index, a measure of market concentration based on
#'     the squared sum of market shares. If \code{scaled = TRUE}, the HHI is normalized to account for the number of firms.}
#'   \item{EG}{Calculates the Ellison-Glaeser Index, which measures the geographic concentration of an industry
#'     while controlling for firm size distribution and random distribution effects. It requires data on regions,
#'     industries, and an indicator (e.g., employment or output).}
#' }
#'
#' @return
#' \describe{
#'   \item{LQ}{A data frame with the region column and calculated location quotients for the specified columns,
#'     optionally grouped by \code{.by}.}
#'   \item{HHI}{A numeric value representing the HHI, either scaled or unscaled.}
#'   \item{EG}{A tibble with two columns: the industry name and the corresponding EG index value.}
#' }
#'
#' @examples
#' # Example data
#' data = data.frame(
#'   region = c("National", "Region_A", "Region_B"),
#'   total = c(10000, 4000, 6000),
#'   industry1 = c(2000, 1000, 1000),
#'   industry2 = c(6000, 2000, 4000)
#' )
#'
#' # Calculate Location Quotient
#' LQ(data, region, total, starts_with("industry"))
#'
#' # Calculate HHI
#' x = c(50, 30, 20)
#' # Calculate the raw HHI
#' HHI(x)
#' # Calculate the standard (scaled) HHI
#' HHI(x, scaled = TRUE)
#'
#' # Example data for EG
#' eg_data = data.frame(
#'  region = c("R1", "R1", "R1", "R1", "R2", "R2", "R3", "R3", "R1", "R2", "R2", "R3"),
#'  industry = c("A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B"),
#'  employment = c(250, 200, 150, 100, 20, 15, 10, 5, 50, 200, 150, 50)
#' )
#' EG(eg_data, region, industry, y = employment)
#'
#' @name regional_economics
NULL

#' @rdname regional_economics
#' @export
LQ = function(data, region, total, .cols, .by = NULL) {
  # Calculate location quotient for multiple columns with grouping
  # data: a data frame, where the first row is the national total, the first two columns are the region column and the total column,
  #       optionally includes a grouping column, and the remaining columns are used to calculate the location quotient
  # region, total: specify the region and total columns, respectively
  # .cols: specifies the columns for which to calculate the location quotient
  # .by: specifies the grouping column, defaults to `NULL` (no grouping)
  quotient = function(x) x[-1] / x[1]
  data |>
    dplyr::reframe({{region}} := {{region}}[-1],
                   dplyr::across({{.cols}}, function(x) quotient(x / {{total}})),
                   .by = {{.by}})
}

#' @rdname regional_economics
#' @export
HHI = function(x, scaled = FALSE) {
  # x: a numeric vector
  # scaled: a logical value, indicating whether to standardize the HHI
  x = x[!is.na(x)]
  hhi = sum((x / sum(x)) ^ 2)
  if(scaled) {
    n = length(x)
    hhi = (hhi - 1/n) / (1 - 1/n)
  }
  hhi
}

#' @rdname regional_economics
#' @export
EG = function(data, region, industry, y) {
  # data: a data frame containing region, industry, and indicator columns (e.g., employment or output)
  # region, industry, y: specify the region, industry, and indicator columns, respectively
  cal_eg = function(s, x, h) {
    if(h == 1) 0
    else (sum((s-x)^2) / (1 - sum(x^2)) - h) / (1 - h)
  }
  H = data |>
    dplyr::summarise(h = HHI({{y}}, scaled = FALSE), .by = {{industry}}) |>
    dplyr::arrange({{industry}})
  X = data |>
    dplyr::summarise(x = sum({{y}}), .by = {{region}}) |>
    dplyr::mutate(x = x / sum(x)) |>
    dplyr::arrange({{region}})
  S = data |>
    dplyr::summarise(s = sum({{y}}), .by = c({{region}}, {{industry}})) |>
    dplyr::mutate(s = s / sum(s), .by = {{industry}}) |>
    dplyr::arrange({{industry}}, {{region}}) |>
    tidyr::pivot_wider(names_from = {{industry}}, values_from = s,
                       values_fill = 0)
  purrr::map2_dbl(S[-1], H$h, function(s, h) cal_eg(s, X$x, h)) |>
    tibble::enframe(name = rlang::englue("{{industry}}"), value = "EG")
}
