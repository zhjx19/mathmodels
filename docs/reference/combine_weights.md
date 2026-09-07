# Combine Subjective and Objective Weights

Combines subjective and objective weights using linear, multiplicative,
or game theory-based methods (geometric mean or linear system).

## Usage

``` r
combine_weights(w_subj, w_obj, type = "linear", alpha = 0.5)
```

## Arguments

- w_subj:

  Numeric vector of subjective weights.

- w_obj:

  Numeric vector of objective weights.

- type:

  Character string specifying the combination method: "linear",
  "multiplicative", "game", or "game_linear".

- alpha:

  Numeric value between 0 and 1, used only for the linear method to
  weight subjective weights. Defaults to 0.5.

## Value

A numeric vector of combined weights. For `multiplicative` and `game`,
the result is normalized to sum to 1. For `linear` and `game_linear`,
the theoretically-defined combination is returned as is: it sums to 1
when both input weight vectors do (only the two normalized methods are
rescaled).

## Details

The function supports four methods:

- Linear: Combines weights as alpha \* w_subj + (1 - alpha) \* w_obj.

- Multiplicative: Combines weights as w_subj \* w_obj, requiring
  positive weights, then normalizes to sum to 1.

- Game: Uses the geometric mean (sqrt(w_subj \* w_obj)) to balance
  weights, then normalizes to sum to 1.

- Game_linear: Uses a game-theoretic approach by solving a linear system
  based on the cross-product of weights. Inputs are expected to be
  weight vectors that each sum to 1.

## Examples

``` r
w_subj = c(0.4, 0.3, 0.2, 0.1)
w_obj = c(0.25, 0.2, 0.3, 0.25)
combine_weights(w_subj, w_obj, type = "linear", alpha = 0.6)
#> [1] 0.34 0.26 0.24 0.16
combine_weights(w_subj, w_obj, type = "multiplicative")
#> [1] 0.4081633 0.2448980 0.2448980 0.1020408
combine_weights(w_subj, w_obj, type = "game")
#> [1] 0.3279556 0.2540333 0.2540333 0.1639778
combine_weights(w_subj, w_obj, type = "game_linear")
#> [1] 0.3735683 0.2823789 0.2176211 0.1264317
```
