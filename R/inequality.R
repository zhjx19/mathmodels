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
#' For \code{theil}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Name of population variable (character).
#' @param y For \code{gini0}, \code{gini}, \code{theil0}: Numeric vector of individual incomes.
#' @param pop For \code{gini}: Numeric vector of group populations or population shares.
#' For \code{theil}, \code{theil0_g}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Name of income variable (character).
#' @param data For \code{theil0_g}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: Data frame containing variables.
#' @param group For \code{theil0_g}, \code{theil_g}: Name of grouping variable (e.g., province).
#' @param group1 For \code{theil_g2_cross}, \code{theil_g2_nest}: Name of first grouping variable (e.g., region or province).
#' @param group2 For \code{theil_g2_cross}, \code{theil_g2_nest}: Name of second grouping variable (e.g., type or city).
#'
#' @return
#' For \code{gini0}, \code{gini}: Numeric Gini coefficient (0 to 1).
#' For \code{theil0}, \code{theil}: Numeric Theil index.
#' For \code{theil0_g}, \code{theil_g}, \code{theil_g2_cross}, \code{theil_g2_nest}: List with two vectors:
#' \itemize{
#'   \item \code{theil}: theil (Theil index and its decomposition),
#'   \item \code{ratio}: ratio (contribution rates of each component).
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
#' data2 = data |> dplyr::count(g, y, name = "pop")
#' theil(data2$y, data2$pop)
#'
#' # Theil index with grouping (sample data)
#' theil0_g(data, "g", "y")
#'
#' # Theil index with grouping (average data)
#' theil_g(data2, "g", "y", "pop")
#'
#' # Theil index with two-level cross-grouping
#' data3 = data.frame(
#'   industry = c("A", "A", "A", "A", "B", "B", "B", "B"),
#'   area = c("East", "East", "West", "West", "East", "East", "West", "West"),
#'   province = c("Shanghai", "Beijing", "Sichuan", "Yunnan", "Shanghai", "Beijing", "Sichuan", "Yunnan"),
#'   avg_wage = c(500, 400, 80, 60, 300, 250, 50, 40),
#'   emp_num = c(100, 80, 90, 70, 120, 100, 80, 60)
#' )
#' theil_g2_cross(data3, "industry", "area", "avg_wage", "emp_num")
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
  n = length(x)
  sum((x[2:n] - x[1:(n-1)]) * (y[2:n] + y[1:(n-1)]) / 2)
}

#' @rdname inequality
#' @export
gini0 = function(y) {
  # Computes Gini coefficient for inequality
  # y: numeric vector of non-negative values
  n = length(y)
  d = outer(y, y, FUN = function(a, b) abs(a - b))
  sum(d) / (2 * n^2 * mean(y))
}

#' @rdname inequality
#' @export
gini = function(y, pop) {
  # Computes Gini coefficient for grouped data
  # y: numeric vector of group incomes or income shares
  # pop: numeric vector of group populations or population shares
  ord = order(y / pop)
  y = y[ord]
  pop = pop[ord]
  x = c(0, cumsum(pop) / sum(pop))
  z = c(0, cumsum(y) / sum(y))
  1 - 2 * trapz(x, z)
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
  # Returns a list of two elements:
  #   - theil: Total Theil index (T), between-group (Tb), within-group (Tw), and per-group1 within-group Theil indices
  #   - ratio: corresponding contribution rates

  # Rename columns for consistency
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

  # total theil index
  TT = Tb + Tw

  # within-group contribution rates
  Rwi = dfb$Yr * dfw$Tw / TT

  list(theil = c(T=TT, Tb=Tb, Tw=Tw, setNames(dfw$Tw, dfw$group)),
       ratio = c(T=1, Tb=Tb/TT, Tw=Tw/TT, setNames(Rwi, dfw$group)))
}

#' @rdname inequality
#' @export
theil_g = function(data, group, y, pop) {
  # Computes Theil index and decomposition for grouped average data
  # data: data frame
  # group: name of grouping variable (e.g., province)
  # y: name of average income variable
  # pop: name of population variable
  # Returns a list of two elements:
  #   - theil: Total Theil index (T), between-group (Tb), within-group (Tw), and per-group within-group Theil indices
  #   - ratio: Corresponding contribution rates

  # Rename columns for consistency
  vars = c(group, y, pop)
  data = data |>
    dplyr::select(dplyr::all_of(vars)) |>
    stats::setNames(c("group", "y","pop"))

  dfb = data |>
    dplyr::summarise(y = sum(y*pop)/sum(pop), pop = sum(pop), .by = group) |>
    dplyr::mutate(yp = y*pop, Yr = yp/sum(yp))
  Tb = theil(dfb$y, dfb$pop)

  dfw = data |>
    tidyr::nest(.by = group) |>
    dplyr::mutate(Tw = sapply(data, function(d) theil(d$y, d$pop)))
  Tw = sum(dfb$Yr * dfw$Tw)

  # total theil index
  TT = Tb + Tw

  # within-group contribution rates
  Rwi = dfb$Yr * dfw$Tw / TT

  list(theil = c(T=TT, Tb=Tb, Tw=Tw, setNames(dfw$Tw, dfw$group)),
       ratio = c(T=1, Tb=Tb/TT, Tw=Tw/TT, setNames(Rwi, dfw$group)))
}

#' @rdname inequality
#' @export
theil_g2_cross = function(data, group1, group2, y, pop) {
  # Computes Theil index and decomposition for two-level grouped average data
  # group1: name of first grouping variable (e.g., industry)
  # group2: name of second grouping variable (e.g., area: East, Central, West)
  # y: name of average income variable (e.g., average wage)
  # pop: name of population variable (e.g., employment number)
  # Returns a list with two elements:
  #   - theil: Total Theil index (T), between-group1 (Tb), within-group1 (Tw),
  #            between-group2 (Tw_b), within-group2 (Tw_w), and per-group1 within-group Theil indices
  #   - ratio: corresponding contribution rates

  # Standardize column names
  vars = c(group1, group2, y, pop)
  data = stats::setNames(dplyr::select(data, dplyr::all_of(vars)), c("group1", "group2", "y", "pop"))

  # Compute between-group1 Theil index (by group1)
  dfb = dplyr::mutate(
    dplyr::summarise(data, y = sum(y * pop) / sum(pop), p = sum(pop), .by = group1),
    yp = y * p, Yr = yp / sum(yp))
  Tb = theil(dfb$y, dfb$p)

  # Compute within-group1 Theil index with nested decomposition
  dfw = data |>
    tidyr::nest(.by = group1) |>
    dplyr::mutate(
      # Compute between-group2 Theil index within each group1
      T_g2_b = sapply(data, function(d) {
        d_agg = dplyr::summarise(d, y = sum(y * pop) / sum(pop), pop = sum(pop), .by = group2)
        theil(d_agg$y, d_agg$pop)
      }),
      # Compute within-group2 Theil index for each group2 within group1
      T_g2_w = sapply(data, function(d) {
        d_nested = tidyr::nest(d, .by = group2)
        sum(sapply(d_nested$data, function(dr) {
          if (nrow(dr) > 1) theil(dr$y, dr$pop) * sum(dr$pop) / sum(d$pop)
          else 0
        }), na.rm = TRUE)
      }),
      # Total within-group1 Theil index for each group1
      Tw = T_g2_b + T_g2_w
    )
  Tw = sum(dfb$Yr * dfw$Tw)
  T_g2_b = sum(dfb$Yr * dfw$T_g2_b)
  T_g2_w = sum(dfb$Yr * dfw$T_g2_w)

  # Compute total Theil index as Tb + Tw
  TT = Tb + Tw

  # Compute contribution rates
  Rwi = dfb$Yr * dfw$Tw / TT

  # Return results
  list(
    theil = c(T = TT, Tb = Tb, Tw = Tw, Tw_b = T_g2_b, Tw_w = T_g2_w,
              setNames(dfw$Tw, dfw$group1)),
    ratio = c(T = 1, Tb = Tb/TT, Tw = Tw/TT, Tw_b = T_g2_b/TT, Tw_w = T_g2_w/TT,
              setNames(Rwi, dfw$group1))
  )
}

#' @rdname inequality
#' @export
theil_g2_nest = function(data, group1, group2, y, pop) {
  # Computes Theil index and decomposition for two-level nested grouped data
  # group1: name of first grouping variable (e.g., province)
  # group2: name of second grouping variable (e.g., city)
  # y: name of average income variable (e.g., average wage)
  # pop: name of population variable (e.g., employment number)
  # Returns a list with two elements:
  #   - theil: Total Theil index (T), between-group1 (Tb_g1), between-group2 (Tb_g2),
  #            within-group2 (Tw), and per-group1/group2 Theil indices
  #   - ratio: corresponding contribution rates

  # Standardize column names
  vars = c(group1, group2, y, pop)
  data = stats::setNames(dplyr::select(data, dplyr::all_of(vars)), c("group1", "group2", "y", "pop"))

  # Compute total income for weights
  total_income = sum(data$y * data$pop)

  # Compute between-group1 Theil index (by group1)
  df_g1 = dplyr::summarise(data, y = sum(y * pop) / sum(pop), p = sum(pop), .by = group1) |>
    dplyr::mutate(yp = y * p, Yr = yp / total_income)
  Tb_g1 = theil(df_g1$y, df_g1$p)

  # Compute between-group2 Theil index within each group1
  df_g2 = data |>
    dplyr::summarise(y = sum(y * pop) / sum(pop), p = sum(pop), .by = c(group1, group2)) |>
    tidyr::nest(.by = group1) |>
    dplyr::mutate(
      Tb_g2 = sapply(data, function(d) {
        if (nrow(d) > 1) theil(d$y, d$p) else 0
      }),
      Yr = sapply(data, function(d) sum(d$y * d$p) / total_income)
    )
  Tb_g2 = sum(df_g2$Yr * df_g2$Tb_g2)

  # Compute within-group2 Theil index for each group2
  df_within = data |>
    tidyr::nest(.by = c(group1, group2)) |>
    dplyr::mutate(
      Tw = sapply(data, function(d) {
        if (nrow(d) > 1) theil(d$y, d$pop) else 0
      }),
      y = sapply(data, function(d) sum(d$y * d$pop) / sum(d$pop)),
      p = sapply(data, function(d) sum(d$pop))
    ) |>
    dplyr::mutate(Yr = (y * p) / total_income)
  Tw = sum(df_within$Yr * df_within$Tw, na.rm = TRUE)

  # Compute total Theil index as Tb_g1 + Tb_g2 + Tw
  Tb = Tb_g1 + Tb_g2
  TT = Tb + Tw

  # Compute contribution rates
  # 1. Each group1's between-group2 contribution (weighted)
  R_group1 = (df_g2$Yr * df_g2$Tb_g2) / TT

  # 2. Each group1_group2's within-group2 contribution
  Rwi = df_within$Yr * df_within$Tw / TT

  # 3. Each group1's total contribution (sum of within-group2 contributions within each group1)
  df_within_group1 = df_within |>
    dplyr::group_by(group1) |>
    dplyr::summarise(R_group1_total = sum(Rwi, na.rm = TRUE))
  R_group1_total = setNames(df_within_group1$R_group1_total, df_within_group1$group1)

  # Return results
  list(
    theil = c(T = TT, Tb_g1 = Tb_g1, Tb_g2 = Tb_g2, Tw = Tw,
              setNames(df_g2$Tb_g2, df_g2$group1),
              setNames(df_within$Tw, paste(df_within$group1, df_within$group2, sep = "_"))),
    ratio = c(T = 1, Tb_g1 = Tb_g1/TT, Tb_g2 = Tb_g2/TT, Tw = Tw/TT,
              setNames(R_group1, df_g2$group1),  # Each group1's between-group2 contribution
              setNames(Rwi, paste(df_within$group1, df_within$group2, sep = "_")),  # Each group2's within contribution
              setNames(R_group1_total, names(R_group1_total)))  # Each group1's total within contribution
  )
}
