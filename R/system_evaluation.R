#' @title System Evaluation Functions for Coupling and Obstacle Analysis
#'
#' @description
#' These functions provide tools for system-level evaluation in multi-indicator systems:
#'
#' \itemize{
#'   \item \code{coupling_degree()}: Computes coupling degree, coordination index, and coupling coordination degree for subsystems.
#'   \item \code{obstacle_degree()}: Computes obstacle degrees for secondary indicators to identify key constraints in the system, enabling batch processing with tidyverse for grouping and summarization.
#' }
#'
#' @param data A data frame with normalized scores (usually in [0,1]) as columns, .
#' @param w Optional vector of weights for indicators or subsystems; defaults to equal weights if NULL.
#' @param id_cols Optional character vector of column names in \code{data} to preserve as identifiers (not used in calculations).
#'
#' @return
#' A tibble depending on the function:
#' \describe{
#'   \item{coupling_degree}{A tibble with columns:
#'     \itemize{
#'       \item \code{ID}: Identifier columns specified by \code{id_cols} (if provided).
#'       \item \code{coupling}: Coupling Degree (range 0-1).
#'       \item \code{coord}: Coordination Index (range 0-1).
#'       \item \code{coupling_coord}: Coupling Coordination Degree (range 0-1).
#'     }
#'   }
#'   \item{obstacle_degree}{A tibble with:
#'     \itemize{
#'       \item \code{ID}: Identifier columns specified by \code{id_cols} (if provided).
#'       \item Columns for secondary indicator obstacle degrees (\eqn{O_{ij} = (1 - X_{ij}) \cdot w_{ij}}).
#'     }
#'     Suitable for grouping and summarizing (e.g., with tidyverse) to compute primary indicator obstacle degrees (\( U_i \)).
#'   }
#' }
#'
#' @examples
#' # Sample normalized subsystem scores
#' df = data.frame(
#'   s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
#'   s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573)
#' )
#' # Coupling Degree Analysis
#' coupling_degree(df)        # Equal weights
#' coupling_degree(df, c(0.6, 0.4))
#' # Obstacle Degree Analysis
#' obstacle_degree(df)        # Equal weights
#' obstacle_degree(df, c(0.6, 0.4))
#'
#' @name system_evaluation
NULL

#' @rdname system_evaluation
#' @export
coupling_degree = function(data, w = NULL, id_cols = NULL) {
  # Computes coupling degree(CD), coordination index(CI), and coupling coordination degree(CCD)
  # data: data frame with normalized subsystem scores as columns, Can include ID columns that do not participate in the calculation
  # w: vector of subsystem weights
  # id_cols: Non-subsystem columns to preserve (not used in calculation)
  # Returns: coupling (Coupling Degree), coord (Coordination Index), coupling_coord (Coupling Coordination Degree)

  ID = NULL
  if(!is.null(id_cols)) {
    ID = data[id_cols]
    data = data[setdiff(names(data), id_cols)]
  }

  p = ncol(data)
  if(is.null(w)) w = rep(1/p,p)
  row_prod = apply(data, 1, function(x) prod(pmax(x, 1e-10)))
  row_sum = pmax(rowSums(data), 1e-10)   # Prevent division by 0
  C = p * row_prod^(1/p) / row_sum
  TT = as.matrix(data) %*% w |> as.vector()
  D = sqrt(C * TT)
  tibble::tibble(ID, coupling = C, coord = TT, coupling_coord = D)
}

#' @rdname system_evaluation
#' @export
obstacle_degree = function(data, w = NULL, id_cols = NULL, scaled = FALSE) {
  # Compute obstacle degree for each indicator
  # data: normalized data frame
  # w: weights for indicators, default is equal weights
  # id_cols: Non-subsystem columns to preserve (not used in calculation)
  # scaled: Whether to perform row normalization on the obstacle degree
  # Return: Secondary indicator obstacle degrees

  ID = NULL
  if(!is.null(id_cols)) {
    ID = data[id_cols]
    data = data[setdiff(names(data), id_cols)]
  }

  n = nrow(data)
  m = ncol(data)
  if (is.null(w)) w = rep(1/m, m)
  diff_mat = 1 - as.matrix(data)     # Compute deviation (1 - X_ij)
  O = sweep(diff_mat, 2, w, "*")     # Compute (1 - X_ij) * w_ij
  if(scaled) O = O / rowSums(O)
  tibble::tibble(ID, as.data.frame(O))
}
