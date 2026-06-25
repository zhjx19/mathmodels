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
