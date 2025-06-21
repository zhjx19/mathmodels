#' Data Envelopment Analysis (DEA & SBM)
#'
#' Calculates standard or supper DEA/SBM efficiency scores (e.g., CCR or BCC models)
#' and slacks.
#'
#' @param data A data frame. The first column should contain DMU (Decision Making Unit) names/identifiers.
#'             Subsequent columns are input/output variables.
#' @param inputs A numeric vector of column indices or a character vector of column names
#'               indicating input variables.
#' @param outputs A numeric vector of column indices or a character vector of column names
#'                indicating (desirable) output variables.
#' @param ud_outputs Optional. A numeric vector of column indices or a character vector of
#'                   column names indicating undesirable output variables. Defaults to `NULL`.
#' @param orientation Character string. Model orientation: `"io"` for input-oriented (default),
#'                    or `"oo"` for output-oriented.
#' @param rts Character string. Returns to scale assumption: `"vrs"` for variable returns to
#'            scale (default), or `"crs"` for constant returns to scale.
#'
#' @return A list containing two elements:
#'   \item{efficiency}{A numeric vector of efficiency scores for each DMU.}
#'   \item{slack}{A data frame or matrix containing the slack values for inputs/outputs.}
#'
#' @details
#' This function serves as a wrapper around specific models from the `deaR` package.
#'
#' @importFrom deaR make_deadata model_basic model_sbmeff model_supereff model_sbmsupereff efficiencies slacks
#'
#' @examples
#' # Sample data
#' data = data.frame(
#'   DMU = paste0("DMU", 1:5),
#'   input1 = c(10, 20, 15, 25, 30),
#'   input2 = c(5, 8, 7, 10, 12),
#'   output = c(100, 150, 120, 180, 200),
#'   ud_output = c(10, 15, 12, 20, 25)
#' )
#'
#' # Standard DEA
#' result = basic_DEA(data, inputs = 2:3, outputs = 4)
#' result$efficiency
#'
#' # DEA with undesirable outputs
#' result = basic_DEA(data, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiency
#'
#' # Super-efficiency DEA
#' result = super_DEA(data, inputs = 2:3, outputs = 4)
#' result$efficiency
#'
#' # Standard SBM
#' result = basic_SBM(data, inputs = 2:3, outputs = 4)
#' result$efficiency
#'
#' # SBM with undesirable outputs
#' result = basic_SBM(data, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiency
#'
#' # Super-efficiency SBM
#' result = super_SBM(data, inputs = 2:3, outputs = 4)
#' result$efficiency
#'
#' # Note: According to deaR, the super-efficiency model
#' # does not take into account undesirable inputs/outputs.
#'
#' @name DEA
NULL

#' @rdname DEA
#' @export
basic_DEA = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  # Data frame, first column is DMU names
  # Inputs/outputs/ud_outputs: column indices or names
  # Orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # RTS: "vrs" (variable, default), "crs" (constant)
  # Returns efficiency and slack values
  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs,
                           ud_outputs = ud_outputs)
  model = deaR::model_basic(dat, orientation = orientation, rts = rts)
  list(efficiency = deaR::efficiencies(model), slack = deaR::slacks(model))
}

#' @rdname DEA
#' @export
super_DEA = function(data, inputs, outputs,
                     orientation = "io", rts = "vrs") {
  # Data frame, first column is DMU names
  # Inputs/outputs: column indices or names
  # Orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # RTS: "vrs" (variable, default), "crs" (constant)
  # Returns super-efficiency and slack values
  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs)
  model = deaR::model_supereff(dat, orientation = orientation, rts = rts)
  list(efficiency = deaR::efficiencies(model), slack = deaR::slacks(model))
}

#' @rdname DEA
#' @export
basic_SBM = function(data, inputs, outputs, ud_outputs = NULL,
                     orientation = "io", rts = "vrs") {
  # Data frame, first column is DMU names
  # Inputs/outputs/ud_outputs: column indices or names
  # Orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # RTS: "vrs" (variable, default), "crs" (constant)
  # Returns efficiency and slack values
  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs,
                           ud_outputs = ud_outputs)
  model = deaR::model_sbmeff(dat, orientation = orientation, rts = rts)
  list(efficiency = deaR::efficiencies(model), slack = deaR::slacks(model))
}

#' @rdname DEA
#' @export
super_SBM = function(data, inputs, outputs,
                     orientation = "io", rts = "vrs") {
  # Data frame, first column is DMU names
  # Inputs/outputs: column indices or names
  # Orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # RTS: "vrs" (variable, default), "crs" (constant)
  # Returns super-efficiency and slack values
  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs)
  model = deaR::model_sbmsupereff(dat, orientation = orientation, rts = rts)
  list(efficiency = deaR::efficiencies(model), slack = deaR::slacks(model))
}
