#' @title Combine Subjective and Objective Weights
#' @description Combines subjective and objective weights using linear, multiplicative, or game theory-based methods (geometric mean or linear system).
#' @param w_subj Numeric vector of subjective weights.
#' @param w_obj Numeric vector of objective weights.
#' @param type Character string specifying the combination method: "linear", "multiplicative", "game", or "game_linear".
#' @param alpha Numeric value between 0 and 1, used only for the linear method to weight subjective weights. Defaults to 0.5.
#' @return A numeric vector of combined weights, normalized to sum to 1.
#' @details
#' The function supports four methods:
#' - Linear: Combines weights as alpha * w_subj + (1 - alpha) * w_obj.
#' - Multiplicative: Combines weights as w_subj * w_obj, requiring positive weights.
#' - Game: Uses the geometric mean (sqrt(w_subj * w_obj)) to balance weights.
#' - Game_linear: Uses a game-theoretic approach by solving a linear system based on the cross-product of weights.
#' @examples
#' w_subj = c(0.4, 0.3, 0.2, 0.1)
#' w_obj = c(0.25, 0.2, 0.3, 0.25)
#' combine_weights(w_subj, w_obj, type = "linear", alpha = 0.6)
#' combine_weights(w_subj, w_obj, type = "multiplicative")
#' combine_weights(w_subj, w_obj, type = "game")
#' combine_weights(w_subj, w_obj, type = "game_linear")
#'
#' @export
combine_weights = function(w_subj, w_obj, type = "linear", alpha = 0.5) {
  if (any(w_subj < 0) || any(w_obj < 0)) stop("Weights must be non-negative.")

  # Combine weights based on type
  w = switch(type,
             linear = {
               if (alpha < 0 || alpha > 1) stop("Alpha must be between 0 and 1.")
               alpha * w_subj + (1 - alpha) * w_obj
             },
             multiplicative = {
               if (any(w_subj <= 0) || any(w_obj <= 0)) stop("Weights must be positive for multiplicative method.")
               w_subj * w_obj
             },
             game = sqrt(w_subj * w_obj),
             game_linear = {
               # Solve linear system for game theory-based weights
               w_mat = cbind(w_subj, w_obj)
               S = crossprod(w_mat)
               alpha = solve(S, diag(S))
               alpha = alpha / sum(alpha)
               as.vector(w_mat %*% alpha)
             }
  )
  # If needed, perform normalization
  if (type %in% c("multiplicative", "game")) {
    w = w / sum(w)
  }
  w
}
