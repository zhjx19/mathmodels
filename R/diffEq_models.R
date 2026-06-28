# =============================================================================
# mathmodels — ODE solver and classic epidemic / ecological models
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helper: validate and reorder init against expected compartment order
# -----------------------------------------------------------------------------
.check_init = function(init, expected, defaults = NULL) {
  if (is.null(names(init)))
    stop("'init' must be a named numeric vector.")
  if (!is.null(defaults)) {
    for (nm in names(defaults)) {
      if (!nm %in% names(init)) init[nm] = defaults[[nm]]
    }
  }
  missing_vars = setdiff(expected, names(init))
  if (length(missing_vars) > 0)
    stop("'init' is missing compartment(s): ",
         paste(missing_vars, collapse = ", "), ".")
  init[expected]
}

# -----------------------------------------------------------------------------
# Internal helper: build the evaluation environment once and return a closure
# that updates state variables efficiently on each ODE call.
# -----------------------------------------------------------------------------
.build_model_func = function(vars, parsed_eqs, parms) {
  n_vars   = length(vars)
  eq_names = names(parsed_eqs)

  # Pre-allocate a single evaluation environment; parameters are written once.
  eval_env = new.env(parent = baseenv())
  eval_env$t = 0
  for (p in names(parms)) eval_env[[p]] = parms[[p]]
  for (v in vars)          eval_env[[v]] = 0

  function(t, y, parms) {
    # Update time and state variables in place — no new environment allocated.
    eval_env$t = t
    for (i in seq_along(vars)) eval_env[[vars[i]]] = y[[i]]

    dy = vapply(seq_len(n_vars), function(i) {
      tryCatch(
        eval(parsed_eqs[[i]], envir = eval_env),
        error = function(e)
          stop("Error evaluating equation for '", eq_names[i],
               "': ", e$message, call. = FALSE)
      )
    }, numeric(1L))

    list(setNames(dy, vars))
  }
}


# =============================================================================
#' General ODE Solver
#'
#' Solves a system of ordinary differential equations (ODEs) using
#' string-formula equations.  Equation strings are parsed once at call time
#' and reused across all integration steps, so the solver is substantially
#' faster than approaches that call \code{parse()} inside the derivative
#' function.
#'
#' @param y0        A named numeric vector of initial values for all state
#'   variables.
#' @param times     A numeric vector of output times.
#' @param equations A named character vector; names are variable names, values
#'   are derivative expressions (e.g. \code{c(S = "-beta*S*I")}).
#' @param params    A named numeric vector or list of parameters referenced
#'   inside \code{equations}.
#' @param method    Integration method passed to \code{deSolve::ode}.
#'   Default \code{"lsoda"}.
#' @param ...       Additional arguments forwarded to \code{deSolve::ode}.
#'
#' @return A \code{data.frame} with a \code{time} column followed by one column
#'   per state variable.
#'
#' @details
#' All equations are parsed with \code{parse()} exactly once before the
#' integrator starts.  A single pre-allocated environment is reused on every
#' derivative evaluation, avoiding repeated memory allocation and garbage
#' collection.
#'
#' @examples
#' # Built-in wrapper (one-liner)
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' head(sir)
#'
#' # Custom SEIR model with demography
#' seir_demo = ode_solver(
#'   y0    = c(S = 1000, E = 1, I = 0, R = 0),
#'   times = seq(0, 200, by = 1),
#'   equations = c(
#'     S = "mu*N - beta*S*I - mu*S",
#'     E = "beta*S*I - alpha*E - mu*E",
#'     I = "alpha*E - gamma*I - mu*I",
#'     R = "gamma*I - mu*R"
#'   ),
#'   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1, mu = 0.01, N = 1000)
#' )
#' tail(seir_demo)
#'
#' @export
#' @importFrom deSolve ode
ode_solver = function(y0, times, equations, params = NULL,
                       method = "lsoda", ...) {

  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  if (is.null(names(y0)))
    stop("'y0' must be a named numeric vector.", call. = FALSE)
  if (is.null(names(equations)))
    stop("'equations' must be a named character vector.", call. = FALSE)
  if (length(equations) != length(y0))
    stop("'equations' must have the same length as 'y0'.", call. = FALSE)

  extra_vars = setdiff(names(equations), names(y0))
  if (length(extra_vars) > 0)
    stop("'equations' contains variable(s) not found in 'y0': ",
         paste(extra_vars, collapse = ", "), ".", call. = FALSE)

  # ------------------------------------------------------------------
  # 2. Normalise params (before model_func is built)
  # ------------------------------------------------------------------
  if (is.null(params)) {
    params = list()
  } else if (!is.list(params)) {
    if (is.null(names(params)) || any(nchar(names(params)) == 0))
      warning("All parameters should be named for access inside equations.")
    params = as.list(params)
  }

  # ------------------------------------------------------------------
  # 3. Pre-parse all equation strings (done once, not per ODE call)
  # ------------------------------------------------------------------
  vars       = names(y0)
  parsed_eqs = lapply(equations, function(eq) parse(text = eq)[[1]])

  # ------------------------------------------------------------------
  # 4. Build derivative function with pre-allocated evaluation environment
  # ------------------------------------------------------------------
  model_func = .build_model_func(vars, parsed_eqs, params)

  # ------------------------------------------------------------------
  # 5. Integrate
  # ------------------------------------------------------------------
  res = deSolve::ode(
    y      = y0,
    times  = times,
    func   = model_func,
    parms  = params,
    method = method,
    ...
  )

  df = as.data.frame(res)
  colnames(df)[1] = "time"
  df
}


# =============================================================================
#' Malthusian (Exponential) Growth Model
#'
#' Solves the Malthusian growth equation:
#' \deqn{\frac{dN}{dt} = r \, N}
#'
#' The analytical solution is \eqn{N(t) = N_0 \, e^{r t}}, which can be used
#' to verify numerical accuracy.
#'
#' @param init   Named numeric vector, e.g. \code{c(N = 100)}.
#' @param params Named numeric vector, e.g. \code{c(r = 0.3)}.
#'   \describe{
#'     \item{r}{Intrinsic growth rate (per time unit).}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time} and \code{N}.
#'
#' @examples
#' res = model_malthus(
#'   init   = c(N = 100),
#'   params = c(r = 0.3),
#'   times  = seq(0, 10, by = 0.1)
#' )
#' # Compare with analytical solution
#' res$N_exact = 100 * exp(0.3 * res$time)
#' head(res)
#'
#' @export
model_malthus = function(init, params, times, ...) {
  init = .check_init(init, expected = "N")
  ode_solver(init, times,
             equations = c(N = "r * N"),
             params    = params, ...)
}


# =============================================================================
#' Logistic Population Growth Model
#'
#' Solves the logistic growth equation:
#' \deqn{\frac{dN}{dt} = r \left(1 - \frac{N}{K}\right) N}
#'
#' The analytical solution is
#' \eqn{N(t) = K \,/\, \bigl(1 + (K/N_0 - 1)\,e^{-r t}\bigr)},
#' which can be used to verify numerical accuracy.
#'
#' @param init   Named numeric vector, e.g. \code{c(N = 10)}.
#' @param params Named numeric vector, e.g. \code{c(r = 0.5, K = 100)}.
#'   \describe{
#'     \item{r}{Intrinsic growth rate.}
#'     \item{K}{Carrying capacity.}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time} and \code{N}.
#'
#' @examples
#' res = model_logistic(
#'   init   = c(N = 10),
#'   params = c(r = 0.5, K = 100),
#'   times  = seq(0, 20, by = 0.1)
#' )
#' head(res)
#'
#' @export
model_logistic = function(init, params, times, ...) {
  init = .check_init(init, expected = "N")
  ode_solver(init, times,
             equations = c(N = "r * (1 - N/K) * N"),
             params    = params, ...)
}


# =============================================================================
#' SI Epidemic Model
#'
#' Two-compartment model with no recovery:
#' \deqn{\frac{dS}{dt} = -\beta S I}
#' \deqn{\frac{dI}{dt} =  \beta S I}
#'
#' Population size \eqn{N = S + I} is conserved.
#'
#' @param init   Named numeric vector, e.g. \code{c(S = 990, I = 10)}.
#' @param params Named numeric vector, e.g. \code{c(beta = 0.002)}.
#'   \describe{
#'     \item{beta}{Transmission rate (contacts per individual per time).}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{S}, \code{I}.
#'
#' @examples
#' model_si(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002),
#'   times  = seq(0, 30, by = 0.1)
#' )
#'
#' @export
model_si = function(init, params, times, ...) {
  init = .check_init(init, expected = c("S", "I"))
  ode_solver(init, times,
             equations = c(
               S = "-beta * S * I",
               I = "beta * S * I"
             ),
             params = params, ...)
}


# =============================================================================
#' SIS Epidemic Model
#'
#' Two-compartment model with recovery but no lasting immunity:
#' \deqn{\frac{dS}{dt} = -\beta S I + \gamma I}
#' \deqn{\frac{dI}{dt} =  \beta S I - \gamma I}
#'
#' Population size \eqn{N = S + I} is conserved.
#'
#' @param init   Named numeric vector, e.g. \code{c(S = 990, I = 10)}.
#' @param params Named numeric vector, e.g. \code{c(beta = 0.002, gamma = 0.1)}.
#'   \describe{
#'     \item{beta}{Transmission rate.}
#'     \item{gamma}{Recovery rate (1/gamma = mean infectious period).}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{S}, \code{I}.
#'
#' @examples
#' model_sis(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#'
#' @export
model_sis = function(init, params, times, ...) {
  init = .check_init(init, expected = c("S", "I"))
  ode_solver(init, times,
             equations = c(
               S = "-beta * S * I + gamma * I",
               I = "beta * S * I - gamma * I"
             ),
             params = params, ...)
}


# =============================================================================
#' SIR Epidemic Model
#'
#' Three-compartment model with permanent immunity after recovery:
#' \deqn{\frac{dS}{dt} = -\beta S I}
#' \deqn{\frac{dI}{dt} =  \beta S I - \gamma I}
#' \deqn{\frac{dR}{dt} =  \gamma I}
#'
#' Population size \eqn{N = S + I + R} is conserved.
#' The basic reproduction number is \eqn{R_0 = \beta N / \gamma}.
#'
#' @param init   Named numeric vector.  \code{R} defaults to 0 if omitted,
#'   e.g. \code{c(S = 990, I = 10)}.
#' @param params Named numeric vector, e.g. \code{c(beta = 0.002, gamma = 0.1)}.
#'   \describe{
#'     \item{beta}{Transmission rate.}
#'     \item{gamma}{Recovery rate (1/gamma = mean infectious period).}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{S}, \code{I},
#'   \code{R}.
#'
#' @examples
#' model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#'
#' @export
model_sir = function(init, params, times, ...) {
  init = .check_init(init,
                      expected = c("S", "I", "R"),
                      defaults = list(R = 0)
  )
  ode_solver(init, times,
             equations = c(
               S = "-beta * S * I",
               I = "beta * S * I - gamma * I",
               R = "gamma * I"
             ),
             params = params, ...)
}


# =============================================================================
#' SEIR Epidemic Model
#'
#' Four-compartment model that adds an exposed (latent) class:
#' \deqn{\frac{dS}{dt} = -\beta S I}
#' \deqn{\frac{dE}{dt} =  \beta S I - \alpha E}
#' \deqn{\frac{dI}{dt} =  \alpha E - \gamma I}
#' \deqn{\frac{dR}{dt} =  \gamma I}
#'
#' \itemize{
#'   \item \strong{S} Susceptible — not yet infected.
#'   \item \strong{E} Exposed — infected but not yet infectious
#'         (mean latent period \eqn{1/\alpha}).
#'   \item \strong{I} Infectious — actively spreading the pathogen
#'         (mean infectious period \eqn{1/\gamma}).
#'   \item \strong{R} Recovered — permanently immune.
#' }
#'
#' Population size \eqn{N = S + E + I + R} is conserved.
#'
#' @param init   Named numeric vector.  \code{R} defaults to 0 if omitted,
#'   e.g. \code{c(S = 980, E = 10, I = 10)}.
#' @param params Named numeric vector, e.g.
#'   \code{c(beta = 0.3, alpha = 0.2, gamma = 0.1)}.
#'   \describe{
#'     \item{beta}{Transmission rate.}
#'     \item{alpha}{Rate of progression from exposed to infectious
#'           (1/alpha = mean latent period).}
#'     \item{gamma}{Recovery rate (1/gamma = mean infectious period).}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{S}, \code{E},
#'   \code{I}, \code{R}.
#'
#' @examples
#' model_seir(
#'   init   = c(S = 980, E = 10, I = 10),
#'   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#'
#' @export
model_seir = function(init, params, times, ...) {
  init = .check_init(init,
                      expected = c("S", "E", "I", "R"),
                      defaults = list(R = 0)
  )
  ode_solver(init, times,
             equations = c(
               S = "-beta * S * I",
               E = "beta * S * I - alpha * E",
               I = "alpha * E - gamma * I",
               R = "gamma * I"
             ),
             params = params, ...)
}


# =============================================================================
#' Lotka-Volterra Predator-Prey Model
#'
#' Classic two-species predator-prey system:
#' \deqn{\frac{dx}{dt} = \alpha \, x - \beta \, x y}
#' \deqn{\frac{dy}{dt} = \delta \, x y - \mu \, y}
#'
#' \itemize{
#'   \item \strong{x} Prey population.
#'   \item \strong{y} Predator population.
#' }
#'
#' @param init   Named numeric vector, e.g. \code{c(x = 40, y = 9)}.
#' @param params Named numeric vector with four parameters:
#'   \describe{
#'     \item{alpha}{Prey intrinsic growth rate.}
#'     \item{beta}{Predation rate (prey killed per predator per unit time).}
#'     \item{delta}{Predator reproduction rate per prey consumed.}
#'     \item{mu}{Predator mortality rate.  Named \code{mu} (not \code{gamma})
#'           to avoid confusion with the recovery rate used in epidemic models.}
#'   }
#' @param times  Numeric vector of output times.
#' @param ...    Additional arguments passed to \code{ode_solver}.
#'
#' @return A \code{data.frame} with columns \code{time}, \code{x}, \code{y}.
#'
#' @examples
#' model_lv(
#'   init   = c(x = 40, y = 9),
#'   params = c(alpha = 0.4, beta = 0.04, delta = 0.02, mu = 0.5),
#'   times  = seq(0, 100, by = 0.1)
#' )
#'
#' @export
model_lv = function(init, params, times, ...) {
  init = .check_init(init, expected = c("x", "y"))
  ode_solver(init, times,
             equations = c(
               x = "alpha * x - beta * x * y",
               y = "delta * x * y - mu * y"
             ),
             params = params, ...)
}


# =============================================================================
# Epidemic Visualization Toolkit for mathmodels
# =============================================================================
#
# A collection of ggplot2-based visualizations and key metrics for compartmental
# epidemic models.  All functions accept a data frame returned by any modelling
# function (model_sir, model_seir, or a custom ode_solver call).
#
# Functions:
#   plot_compartments()      Compartment trajectories
#   plot_incidence()         Daily new infections (dI) with peak marker
#   plot_phase_si()          S–I phase portrait
#   plot_Rt_estimate()       Effective reproduction number R_t
#   epi_metrics()            Key scalar metrics (R0, peak, attack rate)
#
# All functions require ggplot2; some also require tidyr (for pivot_longer).

# =============================================================================
# Internal helper: convert wide ODE output to long format

#' @importFrom tidyr pivot_longer

.to_long_states = function(df) {
  tidyr::pivot_longer(
    df,
    cols = !"time",
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
# 2. Daily new infections (dI) with peak marker

#' Plot daily new infections (dI)
#'
#' Computes the daily change in the infectious compartment \eqn{\Delta I}
#' directly via \code{diff()} and plots it as a line chart.  The peak is
#' highlighted with a point marker.
#'
#' @param df  A data frame returned by \code{ode_solver()} or any
#'   built-in model function.  Must contain columns \code{time} and
#'   \code{I}.
#'
#' @return A \link[ggplot2:ggplot]{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
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

  if (!"I" %in% names(df)) stop("Requires column I.", call. = FALSE)

  df = df[order(df$time), , drop = FALSE]

  # dI = diff(I), drop first row (NA from differencing)
  dI = diff(df$I)
  inc = data.frame(
    time = df$time[-1],
    dI   = dI
  )

  # Find peak
  peak = inc[which.max(inc$dI), ]

  ggplot2::ggplot(inc, ggplot2::aes(x = .data$time, y = .data$dI)) +
    ggplot2::geom_line(linewidth = 1, color = "#d62728") +
    ggplot2::geom_point(data = peak,
                        ggplot2::aes(x = .data$time, y = .data$dI),
                        size = 4, color = "blue") +
    ggplot2::labs(
      title = "Daily new infections",
      x = "Time",
      y = expression(Delta * I)
    ) +
    ggplot2::theme_minimal()
}



# =============================================================================
# 3. Phase plot (S vs I)
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
    stop("Requires columns S and I.", call. = FALSE)
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
# 4. Effective reproduction number R_t
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
#' @importFrom rlang .data
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

  if (!"S" %in% names(df)) stop("Column 'S' required.", call. = FALSE)
  params = as.list(params)
  if (is.null(params[["beta"]]) || is.null(params[["gamma"]])) {
    stop("params must contain beta and gamma.", call. = FALSE)
  }

  df = df[order(df$time), , drop = FALSE]

  if (is.null(N)) {
    state_cols = setdiff(names(df), "time")
    N = rowSums(df[, state_cols, drop = FALSE])[1]
  }

  beta  = params[["beta"]]
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

# =============================================================================
# 5. Key epidemic metrics
# =============================================================================
#' Extract key epidemic metrics
#'
#' Returns a named list of 4 core scalar metrics derived from ODE output.
#'
#' @param df    A data frame with columns \code{time}, \code{S}, \code{I}.
#' @param beta  Transmission rate (numeric).
#' @param gamma Recovery / removal rate (numeric).
#' @param N     Total population (numeric). If \code{NULL} (default), computed
#'   as the sum of compartment values at the first time point.
#' @return A named list with components: \code{R0} (basic reproduction number),
#'   \code{peak_infection} (maximum number of infectious individuals),
#'   \code{peak_time} (time at which peak occurs),
#'   \code{attack_rate} (proportion of susceptible that became infected).
#' @param gamma Recovery / removal rate (numeric).#' @param N     Total population (numeric). If \code{NULL} (default), computed
#'   as the sum of compartment values at the first time point.
#' @param gamma Recovery / removal rate (numeric).
#' @param N     Total population (numeric). If \code{NULL} (default), computed
#'   as the sum of compartment values at the first time point.
#'   \code{attack_rate}.
#'
#' @examples
#' sir = model_sir(
#'   init   = c(S = 990, I = 10),
#'   params = c(beta = 0.002, gamma = 0.1),
#'   times  = seq(0, 50, by = 0.1)
#' )
#' epi_metrics(sir, beta = 0.002, gamma = 0.1)
#'
#' @export
epi_metrics = function(df, beta, gamma, N = NULL) {

  if (!"time" %in% names(df)) stop("Column 'time' required.", call. = FALSE)
  if (!"S"    %in% names(df)) stop("Column 'S' required.",    call. = FALSE)
  if (!"I"    %in% names(df)) stop("Column 'I' required.",    call. = FALSE)
  if (missing(beta)  || is.null(beta))  stop("Must provide beta.",  call. = FALSE)
  if (missing(gamma) || is.null(gamma)) stop("Must provide gamma.", call. = FALSE)

  if (is.null(N)) {
    state_cols = setdiff(names(df), "time")
    N = rowSums(df[, state_cols, drop = FALSE])[1]
  }

  list(
    R0             = beta * N / gamma,
    peak_infection = max(df$I, na.rm = TRUE),
    peak_time      = df$time[which.max(df$I)],
    attack_rate    = (df$S[1] - df$S[nrow(df)]) / N
  )
}
