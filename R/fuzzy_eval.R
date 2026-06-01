#' Fuzzy Comprehensive Evaluation
#'
#' @description
#' Performs fuzzy comprehensive evaluation using different fuzzy composition operators
#' to combine factor weights with a fuzzy evaluation matrix. Suitable for multi-criteria
#' decision analysis with weights from methods like AHP, entropy, CRITIC, CV, or PCA.
#'
#' @param w Numeric vector, factor weights (e.g., from `combine_weights_linear`).
#' @param R Numeric matrix, fuzzy evaluation matrix with columns as factors and rows
#' as evaluation grades. Values should be in [0, 1].
#' @param type Integer or character (1-5), specifying the fuzzy composition operator:
#' \itemize{
#'   \item 1: Min-max (main factor decisive).
#'   \item 2: Product-max (main factor prominent).
#'   \item 3: Weighted sum (additive average).
#'   \item 4: Bounded sum of mins (min-sum bounded).
#'   \item 5: Normalized min-sum (balanced average).
#' }
#'
#' @return A numeric vector of normalized comprehensive evaluation results, summing to 1.
#'
#' @details
#' The function computes a fuzzy comprehensive evaluation vector `B` based on
#' the weight vector `w` and fuzzy evaluation matrix `R`. Five composition
#' operators are supported:
#' \itemize{
#'   \item Type 1 (min-max): `max(min(w, R[j,]))`, emphasizes the main factor.
#'   \item Type 2 (product-max): `max(w * R[j,])`, highlights the main factor.
#'   \item Type 3 (weighted sum): `sum(w * R[j,])`, additive average.
#'   \item Type 4 (bounded sum): `min(1, sum(min(w, R[j,])))`, bounds the sum of mins.
#'   \item Type 5 (normalized min-sum): `sum(min(w, R[j,]/sum(R[j,])))`, balanced average.
#' }
#' The output `B` is normalized to sum to 1. If the sum is zero, an error is thrown.
#' Uses base R for lightweight implementation.
#'
#' @examples
#' w = c(0.3, 0.3, 0.3, 0.1)  # weights (e.g., from AHP or entropy)
#'
#' # fuzzy evaluation matrix (3 grades for 4 factors)
#' R = matrix(c(0.8, 0.7, 0.6, 0.7,
#'              0.1, 0.2, 0.2, 0.1,
#'              0.1, 0.1, 0.2, 0.2), nrow = 3, byrow = TRUE)
#' # Apply fuzzy comprehensive evaluation
#' fuzzy_eval(w, R, type = 3)  # Weighted sum
#'
#' @export
fuzzy_eval = function(w, R, type) {
  # Perform fuzzy comprehensive evaluation
  # w: factor weights, R: fuzzy evaluation matrix, type: operator type (1-5)

  if(!is.numeric(w) || is.matrix(w))
    stop("w must be a numeric vector.")
  if(!is.matrix(R))
    stop("R must be a matrix.")
  if(length(w) != ncol(R))
    stop("Length of w must equal ncol(R).")
  type = as.character(type)
  if(!type %in% as.character(1:5))
    stop("type must be 1, 2, 3, 4, or 5.")

  m = nrow(R)
  n = ncol(R)
  B = rep(0, m)
  for(j in 1:m) {
    switch(type,
           # Min-max, main factor decisive
           "1" = {
             B[j] = max(apply(rbind(w, R[j,]), 2, min))
           },
           #  Product-max, main factor prominent
           "2" = {
             B[j] = max(w * R[j,])
           },
           # Weighted sum, additive average
           "3" = {
             B[j] = sum(w * R[j,])
           },
           # Bounded sum of mins
           "4" = {
             B[j] = min(c(1, sum(apply(rbind(w, R[j,]), 2, min))))
           },
           # Normalized min-sum, balanced average
           "5" = {
             r0 = sum(R[j,])
             B[j] = sum(apply(rbind(w, R[j,] / r0), 2, min))
           }
    )
  }
  if(sum(B) == 0) stop("Normalization failed: sum of B is zero.")
  B / sum(B)      # Normalize
}
