#' Inequality Indices
#'
#' @description
#' Computes inequality indices:
#' \code{gini0} calculates the Gini coefficient for individual sample data.
#' \code{gini} calculates the Gini coefficient for grouped data using income and population shares.
#' \code{theil0} calculates the Theil index for individual sample data.
#' \code{theil} calculates the Theil index for grouped average data.
#' \code{theil0_g} calculates the Theil index and decomposition for grouped sample data.
#' \code{theil_g} calculates the Theil index and decomposition for grouped average data.
#' \code{theil_g2} calculates the Theil index and decomposition for two-level grouped average data.
#'
#' @param x For \code{gini0}, \code{gini}: Numeric vector of non-negative values (e.g., income).
#' @param pop For \code{gini}: Numeric vector of group populations or population shares.
#' For \code{theil}, \code{theil_g}: Name of population variable (character).
#' For \code{theil_g2}: Name of population variable (character, aliased as `pop`).
#' @param y For \code{theil0}: Numeric vector of individual incomes.
#' For \code{theil}, \code{theil0_g}, \code{theil_g}, \code{theil_g2}: Name of income variable (character).
#' @param data For \code{theil0_g}, \code{theil_g}, \code{theil_g2}: Data frame containing variables.
#' @param group For \code{theil0_g}, \code{theil_g}: Name of grouping variable (e.g., province).
#' @param group1 For \code{theil_g2}: Name of first grouping variable (e.g., region).
#' @param group2 For \code{theil_g2}: Name of second grouping variable (e.g., type).
#'
#' @return
#' For \code{gini0}, \code{gini}: Numeric Gini coefficient (0 to 1).
#' For \code{theil0}, \code{theil}: Numeric Theil index.
#' For \code{theil0_g}, \code{theil_g}: List with total Theil index (\code{theil}), between-group (\code{Tb}), within-group (\code{Tw}), within-group components (\code{Twi}), and contribution rates (\code{Rb}, \code{Rw}, \code{Rwi}).
#' For \code{theil_g2}: List with total Theil index and decomposition (\code{Theil}) and within-group components (\code{Within}).
#'
#' @examples
#' # Sample data
#' income = c(10, 20, 30, 40, 100)
#' pop = c(100, 150, 200, 250, 300)
#'
#' # Gini coefficient (individual data)
#' gini0(income)
#'
#' # Gini coefficient (grouped data)
#' gini(income, pop)
#'
#' data = data.frame(g = c("A","A",rep("B",10),rep("A",6)),
#'                   y = c(10,10,rep(8,4),rep(6,6),rep(4,4),2,2))
#'
#' data2 = data |>
#'   dplyr::count(g, y, name = "p")
#'
#' # Theil index (individual sample)
#' theil0(data$y)
#'
#' # Theil index (grouped average)
#' theil(data2$y, data2$p)
#'
#' # Theil index with grouping (sample data)
#' theil0_g(data, "g", "y")
#'
#' # Theil index with grouping (average data)
#' theil_g(data2, "g", "y", "p")
#'
#' # Theil index with two-level grouping
#' data3 = data.frame(
#'   region = c("Eastern", "Eastern", "Central", "Central", "Western", "Western", "Northeast", "Northeast"),
#'   type = c("Urban", "Rural", "Urban", "Rural", "Urban", "Rural", "Urban", "Rural"),
#'   pop = c(24491, 21854, 12850, 22321, 12423, 23522, 5930, 4823),
#'   per_income = c(13375, 4720, 8809, 2957, 8783, 2379, 8730, 3379)
#' )
#' theil_g2(data3, "region", "type", "per_income", "pop")
#' theil_g2(data3, "type", "region", "per_income", "pop")
#'
#' @importFrom dplyr select summarise mutate group_by distinct pull n bind_rows
#' @importFrom tidyr nest
#' @name inequality
NULL

trapz = function(x, y) {
  # Computes numerical integral using trapezoidal rule
  S = 0
  for(i in 2:length(x)) {
    S = S + (x[i] - x[i-1]) * (y[i-1] + y[i]) / 2
  }
  S
}

#' @rdname inequality
#' @export
gini0 = function(x) {
  # Computes Gini coefficient for inequality
  # x: numeric vector of non-negative values
  n = length(x)
  d = outer(x, x, FUN = \(a, b) abs(a - b))
  sum(d) / (2 * n^2 * mean(x))
}

#' @rdname inequality
#' @export
gini = function(x, pop) {
  # Computes Gini coefficient for grouped data
  # income: numeric vector of group incomes or income shares
  # pop: numeric vector of group populations or population shares
  y = c(0, cumsum(pop) / sum(pop))
  z = c(0, cumsum(x) / sum(x))
  (0.5 - trapz(y, z)) / 0.5
}

#' @rdname inequality
#' @export
theil0 = function(y) {
  # Computes Theil index for individual sample data
  # y: numeric vector of incomes
  Yr = y / mean(y)
  sum(y / sum(y) * ifelse(Yr>0, log(Yr), 0))
}

#' @rdname inequality
#' @export
theil = function(y, pop) {
  # Computes Theil index for grouped average data
  # y: vector of group average incomes
  # pop: vector of group population sizes
  Yb = sum(y * pop) / sum(pop)
  Yr = y / Yb
  sum(pop / sum(pop) * Yr * ifelse(Yr > 0, log(Yr), 0))
}

#' @rdname inequality
#' @export
theil0_g = function(data, group, y) {
  # Computes Theil index and decomposition for grouped sample data
  # data: data frame
  # group: name of grouping variable (e.g., province)
  # y: name of income variable
  # Returns Theil index, between/within components, and contribution rates
  vars = c(group, y)
  data = data |>
    dplyr::select(dplyr::all_of(vars)) |>
    stats::setNames(c("group", "y"))
  dfb = data |>
    dplyr::summarise(y = mean(y), p = dplyr::n(), .by = group) |>
    dplyr::mutate(yp = y * p, Yr = yp / sum(yp))
  Tb = theil(dfb$y, dfb$p)
  dfw = data |>
    tidyr::nest(.by = group) |>
    dplyr::mutate(Tw = sapply(data, \(d) theil0(d$y)))
  Tw = sum(dfb$Yr * dfw$Tw)
  theil = Tb + Tw
  list(theil=theil, Tb=Tb, Tw=Tw, Twi=dfw$Tw, Rb=Tb/theil,
       Rw=Tw/theil, Rwi=dfb$Yr * dfw$Tw / theil)
}

#' @rdname inequality
#' @export
theil_g = function(data, group, y, p) {
  # Computes Theil index and decomposition for grouped average data
  # data: data frame
  # group: name of grouping variable (e.g., province)
  # y: name of average income variable
  # p: name of population variable
  # Returns Theil index, between/within components, and contribution rates
  vars = c(group, y, p)
  data = data |>
    dplyr::select(dplyr::all_of(vars)) |>
    stats::setNames(c("group", "y","p"))
  dfb = data |>
    dplyr::summarise(y = sum(y * p) / sum(p), p = sum(p), .by = group) |>
    dplyr::mutate(yp = y * p, Yr = yp / sum(yp))
  Tb = theil(dfb$y, dfb$p)
  dfw = data |>
    tidyr::nest(.by = group) |>
    dplyr::mutate(Tw = sapply(data, \(d) theil(d$y, d$p)))
  Tw = sum(dfb$Yr * dfw$Tw)
  theil = Tb + Tw
  list(theil=theil, Tb=Tb, Tw=Tw, Twi=dfw$Tw, Rb=Tb/theil,
       Rw=Tw/theil, Rwi=dfb$Yr * dfw$Tw / theil)
}

#' @rdname inequality
#' @export
theil_g2 = function(data, group1, group2, y, pop) {
  # Computes Theil index and decomposition for two-level grouped average data
  # data: data frame
  # group1, group2: names of grouping variables (e.g., region, type)
  # y: name of average income variable
  # pop: name of population variable
  # Returns total Theil index, within-group components, and contribution rates
  vars = c(group1, group2, y, pop)
  data = data |>
    dplyr::select(dplyr::all_of(vars)) |>
    stats::setNames(c("group1","group2","y","pop"))
  # Total Theil index
  Ta = data |>
    dplyr::mutate(Yb = sum(pop * y) / sum(pop),
                  Yr = y / Yb, Nr = pop / sum(pop)) |>
    dplyr::summarise(theil = sum(Nr * Yr * log(Yr))) |>
    dplyr::pull(theil)
  # Within-group Theil indices and contribution rates
  Twi = data |>
    dplyr::mutate(Yib = sum(pop * y) / sum(pop),
                  Nr = pop / sum(pop), Yr = y / Yib,
                  .by = group1) |>
    dplyr::summarise(Tw = sum(Nr * Yr * log(Yr)), Yib = unique(Yib),
                     Ni = sum(pop), .by = group1) |>
    dplyr::mutate(Yb = sum(Yib * Ni) / sum(Ni),
                  Yr = Yib / Yb, Nr = Ni / sum(Ni),
                  D = Nr * Yr * Tw / Ta)
  # Total, within, and between Theil indices and contribution rates
  rlt = Twi |>
    dplyr::summarise(theil = Ta, Tw = sum(Nr * Yr * Tw),
                     Tb = sum(Nr * Yr * log(Yr)))
  Theil = dplyr::bind_rows(theil = rlt, ratio = rlt / Ta, .id = "type")
  # Within-group components
  Within = Twi |>
    dplyr::select(group1, Twi = Tw, Rwi = D)
  list(Theil = Theil, Within = Within)
}
