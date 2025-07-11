#' Data Envelopment Analysis (DEA & SBM)
#'
#' Calculates standard and super-efficiency DEA and SBM models (CCR, BCC, and slacks-based),
#' including efficiency score, slacks, lambdas, targets, returns, and references,
#' with support for undesirable outputs.
#'
#' @param data A data frame. The first column should contain DMU (Decision Making Unit) names/identifiers.
#'             Subsequent columns are input/output variables.
#' @param inputs A numeric vector of column indices or a character vector of column names
#'               indicating input variables.
#' @param outputs A numeric vector of column indices or a character vector of column names
#'                indicating (desirable) output variables.
#' @param ud_outputs Optional. A numeric vector of denoting the position of undesirable outputs in the outputs. Defaults to `NULL`.
#' @param orientation Character string. Model orientation: `"io"` for input-oriented (default),
#'                    or `"oo"` for output-oriented.
#' @param rts Character string. Returns to scale assumption: `"vrs"` for variable returns to
#'            scale (default), or `"crs"` for constant returns to scale.
#'
#' @return A list containing six elements:
#'   \item{efficiencies}{A numeric vector of efficiency scores for each DMU.}
#'   \item{slacks}{A data frame or matrix containing the slack values for inputs/outputs.}
#'   \item{lambdas}{A matrix or data frame containing the intensity variables (weights) for each DMU, representing the contribution of reference DMUs to the efficiency frontier.}
#'   \item{targets}{A data frame or matrix containing the efficient target values for inputs and outputs (including undesirable outputs) for each DMU.}
#'   \item{returns}{A character vector indicating the returns-to-scale (RTS) status of each DMU, such as "crs" (constant), "vrs" (variable), "irs" (increasing), or "drs" (decreasing).}
#'   \item{references}{A matrix or data frame listing the reference DMUs (peers) contributing to the efficiency frontier for each evaluated DMU.}
#'
#' @details
#' This function provides a simplified interface for computing efficiency scores using radial Data
#' Envelopment Analysis (DEA) models (Charnes et al., 1978; Banker et al., 1984) and Slacks-Based
#' Measure (SBM) models (Tone, 2001) via the `deaR` package. It supports both standard and
#' super-efficiency models under constant (CRS) or variable (VRS) returns to scale, with optional
#' handling of undesirable outputs. The package includes four functions: `basic_DEA`, `super_DEA`,
#' `basic_SBM`, and `super_SBM`, each tailored to specific DEA or SBM variants.
#'
#' - **Model Types**:
#'   - `basic_DEA`: Implements standard radial DEA models, including CCR (CRS) and BCC (VRS),
#'     optimizing radial efficiency measures (input contraction or output expansion).
#'   - `super_DEA`: Implements super-efficiency radial DEA, excluding the evaluated DMU from the
#'     reference set to allow efficiency scores beyond 1 (radial, output-oriented) or below 1
#'     (radial, input-oriented) for efficient DMUs (Andersen & Petersen, 1993).
#'   - `basic_SBM`: Implements standard SBM, optimizing input and output slacks directly for a
#'     non-radial efficiency measure.
#'   - `super_SBM`: Implements super-efficiency SBM, combining SBM with super-efficiency properties,
#'     excluding the evaluated DMU from the reference set.
#'
#' - **Orientation**:
#'   - Input-oriented (`"io"`): Minimizes inputs while maintaining output levels. Efficiency scores
#'     are in (0, 1] ($\theta \leq 1$ for radial models, $\rho$ or $\delta \leq 1$ for SBM models).
#'   - Output-oriented (`"oo"`): Maximizes outputs for given inputs. Efficiency scores are in
#'     (0, 1]. Radial models output $\eta \geq 1$, which is converted to $1/\eta$; SBM models output $1/\rho^*$ or
#'     $1/\delta$ directly.
#'
#' - **Undesirable Outputs**:
#'   - Supported in `basic_DEA`, `super_DEA`, and `basic_SBM` using directional distance functions
#'     (DDF) with direction vector (g_y, -g_b) to increase desirable outputs and decrease undesirable
#'     outputs (Färe & Grosskopf, 2004).
#'   - For `super_SBM`, undesirable outputs are supported by adapting the SBM super-efficiency model
#'     with DDF, ensuring undesirable outputs are minimized appropriately.
#'
#' - **NA Values in Super-Efficiency Models**:
#'   - Super-efficiency models (`super_DEA`, `super_SBM`) may return NA for efficient DMUs due to
#'     infeasible linear programming solutions, particularly under VRS (Andersen & Petersen, 1993).
#'     This occurs when the reference set, excluding the evaluated DMU, cannot form a feasible
#'     production possibility set.
#'   - Users can handle NA values externally, as described in the package vignette. Common approaches
#'     include replacing NA with standard efficiency scores (from `basic_DEA` or `basic_SBM`) or
#'     excluding DMUs with NA values from further analysis.
#'
#' - **Returns to Scale (RTS)**:
#'   - CRS (`"crs"`): Assumes constant returns to scale, suitable for long-run analysis where scale
#'     effects are absent.
#'   - VRS (`"vrs"`): Assumes variable returns to scale, allowing increasing ("irs") or decreasing
#'     ("drs") returns, determined by the sum of intensity variables ($\lambda$).
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
#' result$efficiencies
#'
#' # DEA with undesirable outputs
#' result = basic_DEA(data, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiencies
#'
#' # Super-efficiency DEA
#' result = super_DEA(data, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # Super-efficiency DEA with undesirable outputs
#' result = super_DEA(data, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiencies
#'
#' # Standard SBM
#' result = basic_SBM(data, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # SBM with undesirable outputs
#' result = basic_SBM(data, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#' result$efficiencies
#'
#' # Super-efficiency SBM
#' result = super_SBM(data, inputs = 2:3, outputs = 4)
#' result$efficiencies
#'
#' # Note: According to deaR, the SBM super-efficiency model
#' # does not take into account undesirable inputs/outputs.
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
  # Note: NA in super-efficiency indicates an infeasible solution. Consider using the standard efficiency value as an alternative.
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
  # Orientation: "io" (input-oriented, default), "oo" (output-oriented)
  # RTS: "vrs" (variable, default), "crs" (constant)
  # Returns super-efficiency and slack values
  dat = deaR::make_deadata(data, inputs = inputs, outputs = outputs)
  model = deaR::model_sbmsupereff(dat, orientation = orientation, rts = rts)
  result = summary(model, returnList = TRUE)
  result$efficiencies = setNames(result$efficiencies$eff, result$efficiencies$DMU)
  result
}
