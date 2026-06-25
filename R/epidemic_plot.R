# =============================================================================
# Epidemic Visualization Toolkit for mathmodels
# =============================================================================
#
# A collection of ggplot2-based visualizations for compartmental epidemic
# models.  All functions accept a data frame returned by any modelling
# function (model_sir, model_seir, or a custom ode_solver call).
#
# Functions:
#   plot_compartments()      Compartment trajectories
#   compute_incidence()      New infections / new deaths
#   plot_incidence()         Incidence over time
#   plot_infectious_curve()  Infectious population
#   plot_cumulative_infection()  Cumulative attack size
#   plot_phase_si()          S–I phase portrait
#   plot_Rt_estimate()       Effective reproduction number R_t
#
# All functions require ggplot2; some also require tidyr (for pivot_longer).

# =============================================================================
# Internal helper: convert wide ODE output to long format

# =============================================================================
#' @importFrom tidyr pivot_longer
utils::globalVariables(c(".data", "time"))
.to_long_states = function(df) {
  tidyr::pivot_longer(
    df,
    cols = -time,
    names_to = "compartment",
    values_to = "value"
  )
}

# =============================================================================
# 1. Compartment trajectories over time
# =============================================================================
#' Plot compartment trajectories
#'
#' Draws one line per compartment over time.  Use \code{compartments} to
#' select a subset of the state variables; the default is to plot all
#' columns except \code{time}.
#'
#' @param df           A data frame returned by \code{ode_solver()} or any
#'   built-in model function.
#' @param compartments Character vector of compartment names to plot.  When
#'   \code{NULL} (the default) every column except \code{time} is plotted.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_compartments(sir)
#' plot_compartments(sir, compartments = c("S", "I"))
#'
#' @export
plot_compartments = function(df, compartments = NULL) {

  if (is.null(compartments)) {
    compartments = setdiff(names(df), "time")
  }

  df_long = df[, c("time", compartments)]
  df_long = .to_long_states(df_long)

  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = .data$time, y = .data$value, color = .data$compartment)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Compartment trajectories",
      x = "Time",
      y = "Population"
    ) +
    ggplot2::theme_minimal()
}

# =============================================================================
# 2. Incidence (new infections / deaths)
# =============================================================================
#' Compute incidence from ODE output
#'
#' Estimates daily new infections and new deaths from the change in
#' compartment sizes between consecutive time steps.
#'
#' @param df              A data frame with a \code{time} column and
#'   compartment columns.
#' @param infection_col   Name of the infectious compartment column
#'   (default \code{"I"}).
#' @param death_col       Name of the deceased compartment column
#'   (default \code{"D"}).  If the column does not exist, death
#'   incidence is silently skipped.
#' @param population_col  Name of the susceptible compartment column
#'   (default \code{"S"}).  The negative difference of \code{S} is used as
#'   a second estimate of new infections.
#'
#' @return A data frame with columns \code{time}, \code{new_infection}
#'   (from S depletion), \code{new_infection_I} (from I change), and
#'   \code{new_death} (if \code{death_col} exists).  The first row of each
#'   difference column is \code{NA}.
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' inc = compute_incidence(sir)
#' head(inc)
#'
#' @export
compute_incidence = function(df,
                              infection_col = "I",
                              death_col = "D",
                              population_col = "S") {

  df = df[order(df$time), , drop = FALSE]

  out = data.frame(time = df$time)

  # New infections from S depletion
  if (population_col %in% names(df)) {
    out$new_infection = c(NA, -diff(df[[population_col]]))
  }

  # Change in infectious compartment
  if (infection_col %in% names(df)) {
    out$new_infection_I = c(NA, diff(df[[infection_col]]))
  }

  # Deaths if available
  if (death_col %in% names(df)) {
    out$new_death = c(NA, diff(df[[death_col]]))
  }

  out
}

#' Plot incidence curves
#'
#' Calls \code{compute_incidence()} and draws one line for each
#' incidence metric (new infections from S, new infections from I,
#' new deaths).
#'
#' @param df  A data frame returned by \code{ode_solver()} or any
#'   built-in model function.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_incidence(sir)
#'
#' @export
plot_incidence = function(df) {

  inc = compute_incidence(df)

  inc_long = tidyr::pivot_longer(
    inc,
    cols = !"time",
    names_to = "type",
    values_to = "value"
  )

  ggplot2::ggplot(inc_long,
                  ggplot2::aes(x = .data$time, y = .data$value, color = .data$type)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::labs(
      title = "Incidence over time",
      x = "Time",
      y = "Daily change"
    ) +
    ggplot2::theme_minimal()
}

# =============================================================================
# 3. Infectious curve
# =============================================================================
#' Plot infectious population
#'
#' Draws a single line showing the number of infectious individuals
#' over time.  Requires a column named \code{"I"}.
#'
#' @param df  A data frame with columns \code{time} and \code{I}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_infectious_curve(sir)
#'
#' @export
plot_infectious_curve = function(df) {

  if (!"I" %in% names(df)) {
    stop("Column 'I' not found.")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$I)) +
    ggplot2::geom_line(linewidth = 1, color = "#d62728") +
    ggplot2::labs(
      title = "Infectious population over time",
      x = "Time",
      y = "I"
    ) +
    ggplot2::theme_minimal()
}

# =============================================================================
# 4. Cumulative infection (attack size)
# =============================================================================
#' Plot cumulative infections
#'
#' Computes cumulative infections as \eqn{S(0) - S(t)} and plots the
#' result over time.  Requires a column named \code{"S"}.
#'
#' @param df  A data frame with columns \code{time} and \code{S}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_cumulative_infection(sir)
#'
#' @export
plot_cumulative_infection = function(df) {

  if (!"S" %in% names(df)) {
    stop("Column 'S' not found.")
  }

  df$cum_infection = df$S[1] - df$S

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$cum_infection)) +
    ggplot2::geom_line(linewidth = 1, color = "#1f77b4") +
    ggplot2::labs(
      title = "Cumulative infections",
      x = "Time",
      y = "Cumulative infected"
    ) +
    ggplot2::theme_minimal()
}

# =============================================================================
# 5. Phase plot (S vs I)
# =============================================================================
#' Phase plot S vs I
#'
#' Draws a trajectory in the \eqn{(S, I)} plane.  Useful for
#' visualising the epidemic orbit and threshold behaviour.
#'
#' @param df  A data frame with columns \code{time}, \code{S}, and
#'   \code{I}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_path labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_phase_si(sir)
#'
#' @export
plot_phase_si = function(df) {

  if (!all(c("S", "I") %in% names(df))) {
    stop("Requires columns S and I.")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$S, y = .data$I)) +
    ggplot2::geom_path(linewidth = 1, color = "#2ca02c") +
    ggplot2::labs(
      title = "Phase plot: S vs I",
      x = "Susceptible",
      y = "Infectious"
    ) +
    ggplot2::theme_minimal()
}

# =============================================================================
# 6. Effective reproduction number R_t
# =============================================================================
#' Plot effective reproduction number R_t
#'
#' Estimates the time-varying effective reproduction number
#' \eqn{R_t} directly from compartment data.  Two methods are
#' available:
#' \describe{
#'   \item{\code{"mechanistic"}}{Uses \eqn{R_t = \beta S(t) / \gamma}.
#'     This is the standard formula for an SIR-type model and is
#'     recommended when the model dynamics match the SIR structure.}
#'   \item{\code{"normalized"}}{Uses \eqn{R_t = R_0 \, S(t) / N} with
#'     \eqn{R_0 = \beta N / \gamma}.  Suitable when different
#'     normalisations are desired.}
#' }
#'
#' @param df     A data frame with columns \code{time} and \code{S}.
#' @param params A named numeric vector or list containing at least
#'   \code{beta} and \code{gamma}.
#' @param N      Total population size.  If \code{NULL} (the default),
#'   it is estimated from the sum of all state columns at the first
#'   time step.
#' @param method Estimation method: \code{"mechanistic"} (default) or
#'   \code{"normalized"}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object with a dashed
#'   horizontal line at \eqn{R_t = 1}.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs theme_minimal
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' plot_Rt_estimate(sir, params = c(beta = 0.002, gamma = 0.1))
#'
#' @export
plot_Rt_estimate = function(df,
                             params,
                             N = NULL,
                             method = c("mechanistic", "normalized")) {

  method = match.arg(method)

  if (!"S" %in% names(df)) stop("Column 'S' required.")
  params = as.list(params)
  if (is.null(params[["beta"]]) || is.null(params[["gamma"]])) {
    stop("params must contain beta and gamma.")
  }

  df = df[order(df$time), ]

  if (is.null(N)) {
    state_cols = setdiff(names(df), "time")
    N = rowSums(df[, state_cols, drop = FALSE])[1]
  }

  beta = params[["beta"]]
  gamma = params[["gamma"]]

  if (method == "mechanistic") {
    df$Rt = beta * df$S / gamma
  } else {
    R0 = beta * N / gamma
    df$Rt = R0 * df$S / N
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$Rt)) +
    ggplot2::geom_line(linewidth = 1, color = "#d62728") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::labs(
      title = "Effective reproduction number (R_t)",
      x = "Time",
      y = "R_t"
    ) +
    ggplot2::theme_minimal()
}
