# Markov Chain and Grey-Markov Prediction Models

Implements Markov chain prediction for state sequences and the
Grey-Markov combined prediction model that integrates GM(1,1) with
Markov chain correction.

## Usage

``` r
markov_chain(S, s0, n_steps = 5)

GM11_markov(X, n_ahead = 3, breaks = c(-0.1, -0.05, -0.02, 0, 0.02, 0.05, 0.1))
```

## Arguments

- S:

  For `markov_chain`: State sequence (vector or factor).

- s0:

  For `markov_chain`: Initial state, must be one of the levels in `S`.

- n_steps:

  For `markov_chain`: Number of prediction steps (positive integer,
  default 5).

- X:

  For `GM11_markov`: Numeric vector of original time series data.

- n_ahead:

  For `GM11_markov`: Number of future periods to predict (default 3).

- breaks:

  For `GM11_markov`: Numeric vector of boundaries for state
  classification based on relative error. Note: `-Inf` and `Inf` are
  automatically prepended/appended internally.

## Value

- markov_chain:

  Returns a list with:

  - `trans_mat`: Transition probability matrix.

  - `pred_probs`: Matrix of state probabilities for each step.

  - `pred_states`: Predicted states (most likely) for each step.

  - `pi_final`: Ultimate stationary distribution.

- GM11_markov:

  Returns a data frame with columns:

  - `Period`: Time period labels (T1, T2, ..., T+n).

  - `Raw`: Original values (NA for future periods).

  - `GM11_fitted`: GM(1,1) fitted/predicted values.

  - `err_state`: Relative error state (for historical) or predicted
    state (for future).

  - `adj_eff`: Adjustment effect (midpoint of state interval).

  - `Markov_adj`: Markov chain adjusted values.

## Details

`markov_chain` constructs a transition probability matrix from a state
sequence and performs multi-step prediction. For states that never
appear as a starting state, equal transition probabilities are assigned
to maintain row sums of 1. The ultimate stationary distribution is
computed via eigenvalue decomposition.

`GM11_markov` first fits a GM(1,1) model to the original series, then
classifies relative errors into states using the specified boundaries,
builds a Markov chain model on the error states, and corrects both
historical fitted values and future predictions.

## Examples

``` r
# --- Markov chain prediction ---
# Weather states: Rainy, Cloudy, Sunny
S = factor(c("Sunny", "Sunny", "Cloudy", "Rainy", "Sunny",
             "Cloudy", "Sunny", "Sunny", "Rainy", "Cloudy",
             "Sunny", "Cloudy", "Rainy", "Sunny", "Sunny",
             "Cloudy", "Sunny", "Rainy", "Cloudy", "Sunny"),
            levels = c("Rainy", "Cloudy", "Sunny"))
markov_chain(S, s0 = "Cloudy", n_steps = 3)
#> $trans_mat
#>         
#>              Rainy    Cloudy     Sunny
#>   Rainy  0.0000000 0.5000000 0.5000000
#>   Cloudy 0.3333333 0.0000000 0.6666667
#>   Sunny  0.2222222 0.4444444 0.3333333
#> 
#> $pred_probs
#>        Rainy    Cloudy     Sunny
#> T1 0.3333333 0.0000000 0.6666667
#> T2 0.1481481 0.4629630 0.3888889
#> T3 0.2407407 0.2469136 0.5123457
#> 
#> $pred_states
#> [1] "Sunny"  "Cloudy" "Sunny" 
#> 
#> $pi_final
#>     Rainy    Cloudy     Sunny 
#> 0.2105263 0.3157895 0.4736842 
#> 

# --- Grey-Markov prediction ---
X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
GM11_markov(X, n_ahead = 3)
#> Level ratio test passed!
#>    Period   Raw GM11_fitted err_state adj_eff Markov_adj
#> 1      T1 174.0      174.00         5   0.010     174.00
#> 2      T2 179.0      172.81         6   0.035     174.55
#> 3      T3 183.0      183.94         4  -0.010     182.11
#> 4      T4 189.0      195.78         3  -0.035     189.16
#> 5      T5 207.0      208.38         4  -0.010     206.32
#> 6      T6 234.0      221.80         7   0.075     214.30
#> 7      T7 220.5      236.08         2  -0.075     219.61
#> 8      T8 256.0      251.28         5   0.010     253.82
#> 9      T9 270.0      267.46         5   0.010     270.16
#> 10    T10 285.0      284.68         5   0.010     287.56
#> 11    T+1    NA      303.01         5   0.010     306.07
#> 12    T+2    NA      322.52         5   0.010     325.78
#> 13    T+3    NA      343.29         5   0.010     346.76
```
