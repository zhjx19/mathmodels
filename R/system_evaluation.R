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
#' @param type Either "standard" for the standard coupling formula (results concentrated near 1) or "adjusted" for the revised formula (results more uniformly distributed in [0,1]), as proposed by Wang Shujia, Kong Wei, et al. in "Misconceptions and Corrections of Domestic Coupling Coordination Degree Models, Journal of Natural Resources, 2021, 36(3): 793–810 (In Chinese)"
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
#'   ID = LETTERS[1:6],
#'   s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
#'   s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573),
#'   s3 = c(0.2655, 0.3721, 0.5729, 0.9082, 0.2017, 0.8984)
#' )
#' # Coupling Degree Analysis
#' coupling_degree(df, id_cols = "ID")        # Equal weights
#' coupling_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID",
#'                 type = "adjusted")         # "adjusted" coupling degree
#' # Obstacle Degree Analysis
#' obstacle_degree(df, id_cols = "ID")        # Equal weights
#' obstacle_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID")
#'
#' @name system_evaluation
NULL

#' @rdname system_evaluation
#' @export
coupling_degree = function(data, w = NULL, id_cols = NULL,
                           type = "standard") {
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
  C = switch(type,
    "standard" = {
       row_prod = apply(data, 1, function(x) prod(pmax(x, 1e-10)))
       row_sum = pmax(rowSums(data), 1e-10)   # Prevent division by 0
       p * row_prod^(1/p) / row_sum
     },
     "adjusted" = {
       apply(data, 1, function(x) {
         idx = combn(p, 2)
         C1 = 1 - mean(abs(x[idx[2,]] - x[idx[1,]]))
         C2 = prod(x / max(x))
         sqrt(C1 * C2)
       })
     })
  TT = as.vector(as.matrix(data) %*% w)
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
