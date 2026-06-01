#' Data Envelopment Analysis and Malmquist Productivity Index
#'
#' Computes standard and super-efficiency Data Envelopment Analysis (DEA) models
#' (CCR, BCC, and slacks-based), including efficiency scores, slacks, lambdas,
#' targets, returns to scale, and references, with support for undesirable outputs.
#' Additionally, calculates the Malmquist productivity index to measure productivity
#' changes over time, including efficiency change (EC), technical change (TC), and
#' other decomposition components.
#'
#' @param data For DEA functions (`basic_DEA`, `super_DEA`, `basic_SBM`, `super_SBM`):
#'   A data frame where the first column contains DMU (Decision Making Unit)
#'   names/identifiers, and subsequent columns are input/output variables.
#'   For `malmquist`: A long-format data frame with a period column, DMU identifiers,
#'   and input/output variables.
#' @param inputs A numeric vector of column indices or a character vector of column
#'   names indicating input variables.
#' @param outputs A numeric vector of column indices or a character vector of column
#'   names indicating (desirable) output variables.
#' @param ud_outputs Optional. A numeric vector indicating the position of undesirable
#'   outputs within the `outputs` parameter. Defaults to `NULL`. Not applicable for
#'   `super_SBM` or `malmquist`.
#' @param period For `malmquist` only: A numeric or character index/name indicating
#'   the column in `data` containing time periods.
#' @param orientation Character string. Model orientation: `"io"` (input-oriented,
#'   default for DEA functions) or `"oo"` (output-oriented, default for `malmquist`).
#' @param rts Character string. Returns to scale assumption: `"vrs"` (variable returns
#'   to scale, default) or `"crs"` (constant returns to scale).
#' @param type1 For `malmquist` only: Reference technology for Malmquist index:
#'   `"cont"` (contemporary), `"seq"` (sequential), or `"glob"` (global, default).
#' @param type2 For `malmquist` only: Decomposition method for Malmquist index:
#'   `"fgnz"` (Färe et al., 1994), `"rd"` (Ray and Desli, 1997, default),
#'   `"gl"` (generalized), or `"bias"` (biased).
#'
#' @return
#' For DEA functions (`basic_DEA`, `super_DEA`, `basic_SBM`, `super_SBM`):
#' A list containing six elements:
#' \itemize{
#'   \item `efficiencies`: A named numeric vector of efficiency scores for each
#'     DMU, standardized to (0, 1] for both input- and output-oriented models.
#'   \item `slacks`: A data frame or matrix containing slack values for inputs
#'     and outputs (including undesirable outputs, if specified).
#'   \item `lambdas`: A matrix or data frame of intensity variables
#'     (\eqn{\lambda}{lambda}), representing the contribution of reference DMUs to
#'     the efficiency frontier (self-excluded in super-efficiency models).
#'   \item `targets`: A data frame or matrix of efficient target values for
#'     inputs and outputs, adjusted for undesirable outputs in DDF models.
#'   \item `returns`: A character vector indicating returns-to-scale status for
#'     each DMU: `"crs"` (constant), `"irs"` (increasing), or `"drs"` (decreasing).
#'   \item `references`: A matrix or data frame listing reference DMUs (peers)
#'     contributing to the efficiency frontier (\eqn{\lambda > 0}{lambda > 0}).
#' }
#'
#' For `malmquist`:
#' A data frame containing:
#' \itemize{
#'   \item `Period`: Time period transitions (e.g., "t~t+1").
#'   \item `DMU`: Decision Making Unit identifiers.
#'   \item `mi`: Malmquist productivity index, measuring total productivity change.
#'   \item `ec`: Efficiency change (EC), computed directly or as `mi / tc`
#'     if not available from the model.
#'   \item `tc`: Technical change (TC), measuring frontier shift.
#'   \item `pech`: Pure efficiency change (if applicable, based on decomposition
#'     method).
#'   \item `sech`: Scale efficiency change (if applicable, based on decomposition
#'     method).
#' }
#'
#' @details
#' This package provides a unified interface for computing efficiency scores and
#' productivity changes using the \pkg{deaR} package. It includes five functions:
#' `basic_DEA`, `super_DEA`, `basic_SBM`, `super_SBM`, and
#' `malmquist`, each tailored to specific DEA or productivity analysis models.
#'
#' \itemize{
#'   \item \strong{DEA Models}:
#'     \itemize{
#'       \item `basic_DEA`: Implements standard radial DEA models (CCR for CRS,
#'         BCC for VRS) as described by Charnes et al. (1978) and Banker et al. (1984),
#'         optimizing radial efficiency (input contraction or output expansion).
#'       \item `super_DEA`: Computes super-efficiency radial DEA, excluding the
#'         evaluated DMU from the reference set to allow efficiency scores beyond 1
#'         (output-oriented) or below 1 (input-oriented) for efficient DMUs
#'         (Andersen & Petersen, 1993).
#'       \item `basic_SBM`: Implements standard Slacks-Based Measure (SBM) models
#'         (Tone, 2001), optimizing input and output slacks for a non-radial
#'         efficiency measure.
#'       \item `super_SBM`: Combines SBM with super-efficiency properties,
#'         excluding the evaluated DMU from the reference set. Note: `super_SBM`
#'         does not support undesirable outputs.
#'     }
#'   \item \strong{Malmquist Productivity Index}:
#'     \itemize{
#'       \item `malmquist`: Calculates the Malmquist productivity index to measure
#'         productivity changes over time, decomposing it into efficiency change (EC)
#'         and technical change (TC), with optional further decomposition into pure
#'         efficiency change (PECH) and scale efficiency change (SECH) based on
#'         `type2` (Färe et al., 1994; Ray & Desli, 1997). If EC is unavailable
#'         (e.g., under `rts = "vrs"` and `type1 = "glob"`), it is computed
#'         as `ec = mi / tc` to ensure consistent output.
#'     }
#' }
#'
#' \strong{Orientation}:
#' \itemize{
#'   \item Input-oriented (`"io"`): Minimizes inputs while maintaining outputs.
#'     Efficiency scores are in \eqn{(0, 1]}{(0, 1]} (\eqn{\theta \leq 1}{theta <= 1}
#'     for radial models, \eqn{\rho}{rho} or \eqn{\delta \leq 1}{delta <= 1} for
#'     SBM models).
#'   \item Output-oriented (`"oo"`): Maximizes outputs for given inputs.
#'     Efficiency scores are in \eqn{(0, 1]}{(0, 1]}. Radial models output
#'     \eqn{\eta \geq 1}{eta >= 1}, converted to \eqn{1/\eta}{1/eta}; SBM models
#'     output \eqn{1/\rho^*}{1/rho*} or \eqn{1/\delta}{1/delta}.
#' }
#'
#' \strong{Returns to Scale (RTS)}:
#' \itemize{
#'   \item CRS (`"crs"`): Assumes constant returns to scale, suitable for
#'     long-run analysis.
#'   \item VRS (`"vrs"`): Allows variable returns to scale, with increasing
#'     (`"irs"`) or decreasing (`"drs"`) returns determined by the sum of
#'     intensity variables (\eqn{\lambda}{lambda}).
#' }
#'
#' \strong{Undesirable Outputs}:
#' \itemize{
#'   \item Supported in `basic_DEA`, `super_DEA`, and `basic_SBM`
#'     using directional distance functions (DDF) with direction vector
#'     \eqn{(g_y, -g_b)}{(g_y, -g_b)} to increase desirable outputs and decrease
#'     undesirable outputs (Färe & Grosskopf, 2004). Not supported in `super_SBM`
#'     or `malmquist`.
#' }
#'
#' \strong{Malmquist-Specific Parameters}:
#' \itemize{
#'   \item `type1`: Defines the reference technology for the Malmquist index:
#'     \itemize{
#'       \item `"cont"`: Contemporary technology, using each period's frontier.
#'       \item `"seq"`: Sequential technology, incorporating all prior periods.
#'       \item `"glob"`: Global technology, using a single frontier across all
#'         periods.
#'     }
#'   \item `type2`: Specifies the decomposition method:
#'     \itemize{
#'       \item `"fgnz"`: Färe et al. (1994) decomposition.
#'       \item `"rd"`: Ray and Desli (1997) decomposition.
#'       \item `"gl"`: Generalized decomposition.
#'       \item `"bias"`: Bias-corrected decomposition.
#'     }
#' }
#'
#' \strong{Handling NA Values}:
#' \itemize{
#'   \item Super-efficiency models (`super_DEA`, `super_SBM`) may return
#'     `NA` for efficient DMUs due to infeasible linear programming solutions,
#'     especially under VRS (Andersen & Petersen, 1993). Users can replace `NA`
#'     with standard efficiency scores or exclude affected DMUs.
#'   \item For `malmquist`, if `ec` is `NULL` (e.g., under
#'     `rts = "vrs"` and `type1 = "glob"`), it is computed as
#'     `mi / tc` to ensure a complete result.
#' }
#'
#' The package leverages \pkg{deaR} for robust computation, handling zero values
#' internally and ensuring compatibility with input/output specifications. Efficiency
#' scores are standardized to \eqn{(0, 1]} for DEA models, and Malmquist results are
#' formatted as a data frame for easy analysis.
#'
#' @importFrom deaR make_deadata make_malmquist model_basic model_sbmeff
#' @importFrom deaR model_supereff model_sbmsupereff malmquist_index efficiencies
#' @importFrom deaR slacks
#' @importFrom dplyr mutate
#'
#' @examples
#' # Sample DEA data
#' data_dea = data.frame(
#'   DMU = paste0("DMU", 1:5),
#'   input1 = c(10, 20, 15, 25, 30),
#'   input2 = c(5, 8, 7, 10, 12),
#'   output = c(100, 150, 120, 180, 200),
#'   ud_output = c(10, 15, 12, 20, 25)
#' )
#'
#' # Standard DEA
#' result = basic_DEA(data_dea, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # DEA with undesirable outputs
#' result = basic_DEA(data_dea, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiencies
#'
#' # Super-efficiency DEA
#' result = super_DEA(data_dea, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # Standard SBM
#' result = basic_SBM(data_dea, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # Super-efficiency SBM
#' result = super_SBM(data_dea, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # Sample Malmquist data (long format)
#' data_malm = data.frame(
#'   DMU = rep(paste0("DMU", 1:5), 3),
#'   Period = rep(1:3, each = 5),
#'   input1 = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
#'   input2 = c(5, 8, 7, 10, 12, 6, 9, 8, 11, 13, 7, 10, 9, 12, 14),
#'   output = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
#' )
#' malmquist(data_malm, period = 2, inputs = 3:4, outputs = 5)
#'
#' @name DEA
NULL

#' @rdname DEA
#' @export
basic_DEA = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  # data: Dataframe, first column is DMU names
  # inputs/outputs(contains ud_outputs): column indices or names
  # ud_outputs: integer numbers denoting the position of undesirable outputs in the outputs
  # orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # rts: "vrs" (variable, default), "crs" (constant)
  # returns efficiencies, slacks, lambdas, targets, returns, references

  if(!is.data.frame(data))
    stop("data must be a data frame.")
  if(ncol(data) < 3)
    stop("data must have at least 3 columns (DMU + at least 1 input + at least 1 output).")

  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs,
                           ud_outputs = ud_outputs)
  model = deaR::model_basic(dat, orientation = orientation, rts = rts)
  result = summary(model, returnList = TRUE)
  # Convert efficiency scores for output orientation to (0, 1]
  if (orientation == "oo") {
    result$efficiencies = setNames(1 / result$efficiencies$eff, result$efficiencies$DMU)
  } else {
    result$efficiencies = setNames(result$efficiencies$eff, result$efficiencies$DMU)
  }
  result
}

#' @rdname DEA
#' @export
super_DEA = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  # data: Dataframe, first column is DMU names
  # inputs/outputs(contains ud_outputs): column indices or names
  # ud_outputs: integer numbers denoting the position of undesirable outputs in the outputs
  # orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # rts: "vrs" (variable, default), "crs" (constant)
  # returns efficiencies, slacks, lambdas, targets, returns, references
  # Note: NA in super-efficiency indicates an infeasible solution.

  if(!is.data.frame(data))
    stop("data must be a data frame.")
  if(ncol(data) < 3)
    stop("data must have at least 3 columns (DMU + at least 1 input + at least 1 output).")

  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs,
                           ud_outputs = ud_outputs)
  model = deaR::model_supereff(dat, orientation = orientation, rts = rts)
  result = summary(model, returnList = TRUE)
  # Convert efficiency scores for output orientation to (0, 1]
  if (orientation == "oo") {
    result$efficiencies = setNames(1 / result$efficiencies$eff, result$efficiencies$DMU)
  } else {
    result$efficiencies = setNames(result$efficiencies$eff, result$efficiencies$DMU)
  }
  result
}

#' @rdname DEA
#' @export
basic_SBM = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  # data: Dataframe, first column is DMU names
  # inputs/outputs(contains ud_outputs): column indices or names
  # ud_outputs: integer numbers denoting the position of undesirable outputs in the outputs
  # orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # rts: "vrs" (variable, default), "crs" (constant)
  # returns efficiencies, slacks, lambdas, targets, returns, references

  if(!is.data.frame(data))
    stop("data must be a data frame.")
  if(ncol(data) < 3)
    stop("data must have at least 3 columns (DMU + at least 1 input + at least 1 output).")

  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs,
                           ud_outputs = ud_outputs)
  model = deaR::model_sbmeff(dat, orientation = orientation, rts = rts)
  result = summary(model, returnList = TRUE)
  result$efficiencies = setNames(result$efficiencies$eff, result$efficiencies$DMU)
  result
}

#' @rdname DEA
#' @export
super_SBM = function(data, inputs, outputs,
                     orientation = "io", rts = "vrs") {
  # Note: According to deaR, the SBM super-efficiency model does not take into account undesirable inputs/outputs
  # data: Dataframe, first column is DMU names
  # Inputs/outputs: column indices or names
  # orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # rts: "vrs" (variable, default), "crs" (constant)
  # Returns super-efficiency and slack values

  if(!is.data.frame(data))
    stop("data must be a data frame.")
  if(ncol(data) < 3)
    stop("data must have at least 3 columns (DMU + at least 1 input + at least 1 output).")

  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs)
  model = deaR::model_sbmsupereff(dat, orientation = orientation, rts = rts)
  result = summary(model, returnList = TRUE)
  result$efficiencies = setNames(result$efficiencies$eff, result$efficiencies$DMU)
  result
}

#' @rdname DEA
#' @export
malmquist = function(data, period, inputs, outputs,
                     orientation = "oo", rts = "vrs",
                     type1 = "glob", type2 = "rd") {
  # data: long-format data frame
  # period, inputs, outputs: column indices or names for period, inputs, outputs
  # orientation: "io" (input-oriented), "oo" (output-oriented, default)
  # rts: "vrs" (variable, default), "crs" (constant)
  # type1: reference technology: "cont" (contemporary), "seq" (sequential) or "glob" (global, default)
  # type2: decomposition method: "fgnz" (Fare et al. 1994), "rd" (Ray and Desli 1997, default), "gl" (generalized) or "bias" (biased)
  # Return: data frame with Period, DMU, Malmquist Index, EC, TC, PECH, SECH

  if(!is.data.frame(data))
    stop("data must be a data frame.")
  if(ncol(data) < 4)
    stop("data must have at least 4 columns (period + DMU + at least 1 input + at least 1 output).")

  dat = deaR::make_malmquist(data, percol = period, arrangement = "vertical",
                             inputs = inputs, outputs = outputs)
  model = deaR::malmquist_index(dat, orientation = orientation,
                                rts = rts, type1 = type1, type2 = type2)
  res = summary(model)$Results |>
    dplyr::mutate(Period = paste0(as.numeric(Period) - 1, "~", Period))
  if(is.null(model$ec)) {
    res = res |>
      dplyr::mutate(ec = mi / tc, .after = mi)
  }
  tibble::as_tibble(res)
}
