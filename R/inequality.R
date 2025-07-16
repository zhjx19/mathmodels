#' Inequality Indices
#'
#' @description
#' Computes inequality indices for individual or grouped data:
#' \code{gini0} calculates the Gini coefficient for individual sample data.
#' \code{gini} calculates the Gini coefficient for grouped data using income and population shares.
#' \code{theil0} calculates the Theil index for individual sample data.
#' \code{theil} calculates the Theil index for grouped average data.
#' \code{theil0_g} calculates the Theil index and decomposition for grouped sample data.
#' \code{theil_g} calculates the Theil index and decomposition for grouped average data.
#' \code{theil_g2_cross} calculates the Theil index and decomposition for two-level cross-grouped average data.
#' \code{theil_g2_nest} calculates the Theil index and decomposition for two-level nested grouped average data.
#'
#' @param x For \code{gini0}, \code{gini}: Numeric vector of non-negative values (e.g., income).
#' @param pop For \code{gini}: Numeric vector of group populations or population shares.
#' For \code{theil}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Name of population variable (character).
#' @param y For \code{theil0}: Numeric vector of individual incomes.
#' For \code{theil}, \code{theil0_g}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Name of income variable (character).
#' @param data For \code{theil0_g}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Data frame containing variables.
#' @param group For \code{theil0_g}, \code{theil_g}: Name of grouping variable (e.g., province).
#' @param group1 For \code{theil_g2_cross}, \code{theil_g2_nest}: Name of first grouping variable (e.g., region or province).
#' @param group2 For \code{theil_g2_cross}, \code{theil_g2_nest}: Name of second grouping variable (e.g., type or city).
#'
#' @return
#' For \code{gini0}, \code{gini}: Numeric Gini coefficient (0 to 1).
#' For \code{theil0}, \code{theil}: Numeric Theil index.
#' For \code{theil0_g}, \code{theil_g}: List with two tibbles:
#' \itemize{
#'   \item \code{total}: Tibble with columns \code{type} ("value", "rate"), \code{theil} (total Theil index and 1), \code{Tb} (between-group inequality and contribution rate), \code{Tw} (within-group inequality and contribution rate).
#'   \item \code{within}: Tibble with columns \code{group} (grouping variable), \code{Twi} (within-group Theil indices), \code{Rwi} (within-group contribution rates).
#' }
#' For \code{theil_g2_cross}: List with two tibbles:
#' \itemize{
#'   \item \code{total}: Tibble with columns \code{type} ("value", "rate"), \code{theil} (total Theil index and 1), \code{Tb} (between-group1 inequality and contribution rate), \code{Tw} (within-group1 inequality and contribution rate).
#'   \item \code{within}: Tibble with columns \code{group1} (first grouping variable), \code{Twi} (within-group1 Theil indices), \code{Rwi} (within-group1 contribution rates).
#' }
#' For \code{theil_g2_nest}: Tibble with two rows:
#' \itemize{
#'   \item Row 1 (\code{type = "value"}): Columns \code{theil} (total Theil index), \code{Tw} (within-group2 inequality), \code{Tb} (between-group2 inequality), \code{Tb_group1} (between-group1 inequality), \code{Tb_group2} (within-group1 between-group2 inequality).
#'   \item Row 2 (\code{type = "ratio"}): Columns \code{theil} (1), \code{Tw} (within-group2 contribution rate), \code{Tb} (between-group2 contribution rate), \code{Tb_group1} (between-group1 contribution rate), \code{Tb_group2} (within-group1 between-group2 contribution rate).
#' }
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
#' # Theil index (individual sample)
#' data = data.frame(g = c("A","A",rep("B",10),rep("A",6)),
#'                   y = c(10,10,rep(8,4),rep(6,6),rep(4,4),2,2))
#' theil0(data$y)
#'
#' # Theil index (grouped average)
#' data2 = data |> dplyr::count(g, y, name = "p")
#' theil(data2$y, data2$p)
#'
#' # Theil index with grouping (sample data)
#' theil0_g(data, "g", "y")
#'
#' # Theil index with grouping (average data)
#' theil_g(data2, "g", "y", "p")
#'
#' # Theil index with two-level cross-grouping
#' data3 = data.frame(
#'   region = c("Eastern", "Eastern", "Central", "Central", "Western", "Western", "Northeast", "Northeast"),
#'   type = c("Urban", "Rural", "Urban", "Rural", "Urban", "Rural", "Urban", "Rural"),
#'   pop = c(24491, 21854, 12850, 22321, 12423, 23522, 5930, 4823),
#'   per_income = c(13375, 4720, 8809, 2957, 8783, 2379, 8730, 3379)
#' )
#' theil_g2_cross(data3, "region", "type", "per_income", "pop")
#' theil_g2_cross(data3, "type", "region", "per_income", "pop")
#'
#' # Theil index with two-level nested grouping
#' data4 = data.frame(
#'   province = c("A", "A", "A", "A", "B", "B"),
#'   city = c("A1", "A1", "A2", "A2", "B1", "B1"),
#'   industry = c("Manu", "Serv", "Manu", "Serv", "Manu", "Serv"),
#'   y = c(50000, 45000, 60000, 55000, 70000, 65000),
#'   pop = c(10000, 8000, 15000, 12000, 10000, 8000)
#' )
#' theil_g2_nest(data4, "province", "city", "y", "pop")
#'
#' @importFrom dplyr select summarise mutate n bind_rows
#' @importFrom tidyr nest
#' @importFrom tibble tibble
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
  d = outer(x, x, FUN = function(a, b) abs(a - b))
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
    dplyr::mutate(Tw = sapply(data, function(d) theil0(d$y)))

  Tw = sum(dfb$Yr * dfw$Tw)
  TT = Tb + Tw            # total theil index
  total = tibble(type = c("value", "rate"), theil = c(TT, 1),
                 Tb = c(Tb, Tb/TT), Tw = c(Tw, Tw/TT))
  within = tibble(group = dfw$group, Twi = dfw$Tw,
                  Rwi = dfb$Yr * dfw$Tw / TT)

  list(total = total, within = within)
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
    dplyr::mutate(Tw = sapply(data, function(d) theil(d$y, d$p)))

  Tw = sum(dfb$Yr * dfw$Tw)
  TT = Tb + Tw            # total theil index
  total = tibble(type = c("value", "rate"), theil = c(TT, 1),
                 Tb = c(Tb, Tb/TT), Tw = c(Tw, Tw/TT))
  within = tibble(group = dfw$group, Twi = dfw$Tw,
                  Rwi = dfb$Yr * dfw$Tw / TT)

  list(total = total, within = within)
}

#' @rdname inequality
#' @export
theil_g2_cross = function(data, group1, group2, y, pop) {
  # Computes Theil index and decomposition for two-level grouped data
  # data: data frame
  # group1: name of first grouping variable (e.g., region)
  # group2: name of second grouping variable (e.g., type)
  # y: name of average income variable
  # pop: name of population variable
  # Returns list with Theil index, between/within components, and contribution rates
  # Rename columns for consistency
  vars = c(group1, group2, y, pop)
  data = data |>
    dplyr::select(dplyr::all_of(vars)) |>
    stats::setNames(c("group1", "group2", "y", "pop"))

  # Compute total Theil index using theil()
  Ta = theil(data$y, data$pop)

  # Compute group1-level aggregates for between-group inequality
  dfb = data |>
    dplyr::summarise(y = sum(y * pop) / sum(pop), p = sum(pop),
                     .by = group1) |>
    dplyr::mutate(yp = y * p, Yr = yp / sum(yp))

  # Between-group Theil index using theil()
  Tb = theil(dfb$y, dfb$p)

  # Compute within-group Theil index for each group1
  dfw = data |>
    tidyr::nest(.by = group1) |>
    dplyr::mutate(Tw = sapply(data, function(d) theil(d$y, d$pop)))

  # Total within-group inequality
  Tw = sum(dfb$Yr * dfw$Tw)

  # Contribution rates
  Rwi = dfb$Yr * dfw$Tw / Ta

  # Combine results
  total = bind_rows(
    theil = tibble::tibble(theil = Ta, Tw = Tw, Tb = Tb),
    ratio = tibble::tibble(theil = Ta/Ta, Tw = Tw/Ta, Tb = Tb/Ta),
    .id = "type")
  within = tibble::tibble(group1 = dfw$group1, Twi = dfw$Tw, Rwi = Rwi)

  list(total = total, within = within)
}

#' @rdname inequality
#' @export
theil_g2_nest = function(data, group1, group2, y, pop) {
  # Computes Theil index and decomposition for two-level nested grouped data
  # data: data frame
  # group1: name of first grouping variable (e.g., province)
  # group2: name of second grouping variable (e.g., city)
  # y: name of average income variable
  # pop: name of population variable
  # Returns tibble with Theil index, between/within components, and contribution rates

  # Standardize column names
  data = data |>
    dplyr::select(dplyr::all_of(c(group1, group2, y, pop))) |>
    stats::setNames(c("group1", "group2", "y", "pop"))

  # total theil index
  TT = theil(data$y, data$pop)

  # Aggregate to city level, calculate Tb (inequality between cities)
  df_group2 = data |>
    dplyr::summarise(y = sum(y * pop) / sum(pop),             # City average wage
                     pop = sum(pop), .by = c(group1, group2)) # City total pop
  Tb = theil(df_group2$y, df_group2$pop)

  # Aggregate to provincial level, calculate Tb_group1
  df_prov = data |>
    dplyr::summarise(y = sum(y * pop) / sum(pop),  # Province average wage
                     pop = sum(pop), .by = group1) # Province total pop
  Tb_g1 = theil(df_prov$y, df_prov$pop)

  # Calculate Tb_group2
  Tb_g2 = Tb - Tb_g1

  # Calculate within-group inequality Tw (inequality within cities)
  df_within = data |>
    tidyr::nest(.by = c(group1, group2)) |>
    dplyr::mutate(Tw = sapply(data, function(d) theil(d$y, d$pop)))
  df_g2_weights = df_group2 |>
    dplyr::mutate(Yr = (y * pop) / sum(y * pop))
  Tw = sum(df_g2_weights$Yr * df_within$Tw)

  tibble::tibble(type = c("value", "ratio"), theil = c(TT, 1),
                 Tw = c(Tw, Tw / TT), Tb = c(Tb, Tb / TT),
                 Tb_group1 = c(Tb_g1, Tb_g1 / TT),
                 Tb_group2 = c(Tb_g2, Tb_g2 / TT))
}
