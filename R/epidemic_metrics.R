# =============================================================================
# Epidemiological decision metrics from ODE outputs
# =============================================================================

#' Compute epidemiological decision metrics
#'
#' This function extracts key epidemiological indicators from
#' compartmental model outputs for decision support.
#'
#' @param df Data frame with time series from ODE solver.
#' Must contain at least columns: time, S, I.
#'
#' @param params Named vector or list with:
#' \describe{
#'   \item{beta}{Transmission rate}
#'   \item{gamma}{Recovery rate}
#' }
#'
#' @param N Population size (optional). If NULL, inferred from initial state.
#'
#' @param threshold Healthcare or policy threshold for I(t)
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{Key scalar metrics}
#'   \item{trajectory_metrics}{Time-varying derived indicators}
#' }
#'
#' @export
epidemic_metrics = function(df,
                             params,
                             N = NULL,
                             threshold = NULL) {

  if (!"time" %in% names(df))
    stop("df must contain 'time' column.")

  if (!all(c("S", "I") %in% names(df)))
    stop("df must contain 'S' and 'I' columns.")

  params = as.list(params)
  if (is.null(params$beta) || is.null(params$gamma))
    stop("params must contain beta and gamma.")

  df = df[order(df$time), ]

  # ------------------------------------------------------------
  # Infer population size
  # ------------------------------------------------------------
  if (is.null(N)) {
    state_cols = setdiff(names(df), "time")
    N = rowSums(df[, state_cols, drop = FALSE])[1]
  }

  beta  = params$beta
  gamma = params$gamma

  # ------------------------------------------------------------
  # R0 and Rt
  # ------------------------------------------------------------
  R0 = beta * N / gamma
  df$Rt = R0 * df$S / N

  # ------------------------------------------------------------
  # Incidence
  # ------------------------------------------------------------
  df$incidence = c(NA, -diff(df$S))

  # ------------------------------------------------------------
  # Growth rate of infection
  # ------------------------------------------------------------
  df$growth_rate = c(NA, diff(log(pmax(df$I, 1e-8))))

  # ------------------------------------------------------------
  # Peak infection
  # ------------------------------------------------------------
  peak_I = max(df$I, na.rm = TRUE)
  peak_time = df$time[which.max(df$I)]

  # ------------------------------------------------------------
  # Attack rate (final size)
  # ------------------------------------------------------------
  attack_rate = (df$S[1] - utils::tail(df$S, 1)) / N

  # ------------------------------------------------------------
  # Time above threshold
  # ------------------------------------------------------------
  time_above_threshold = NA
  if (!is.null(threshold)) {
    time_above_threshold = sum(df$I > threshold, na.rm = TRUE)
  }

  # ------------------------------------------------------------
  # Rt < 1 duration (control period)
  # ------------------------------------------------------------
  control_time = sum(df$Rt < 1, na.rm = TRUE)

  # ------------------------------------------------------------
  # Summary output
  # ------------------------------------------------------------
  summary = list(
    R0 = R0,
    peak_infection = peak_I,
    peak_time = peak_time,
    attack_rate = attack_rate,
    final_susceptible = utils::tail(df$S, 1),
    final_infected = utils::tail(df$I, 1),
    control_time_steps = control_time,
    time_above_threshold = time_above_threshold
  )

  return(list(
    summary = summary,
    trajectory = df
  ))
}
