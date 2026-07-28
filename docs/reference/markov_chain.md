# Markov Chain Prediction

Constructs a transition probability matrix from a state sequence and
performs multi-step Markov chain prediction.

## Usage

``` r
markov_chain(S, s0, n_steps = 5)
```

## Arguments

- S:

  State sequence, supplied as a vector or factor.

- s0:

  Initial state, must be one of the state levels in `S`.

- n_steps:

  Number of prediction steps. Must be a positive integer.

## Value

A list with components:

- trans_mat:

  Transition probability matrix.

- pred_probs:

  Matrix of state probabilities for each prediction step.

- pred_states:

  Predicted states, taking the most likely state at each step.

- pi_final:

  Ultimate stationary distribution.

## Details

For states that never appear as a starting state, equal transition
probabilities are assigned to keep each row sum equal to 1. The ultimate
stationary distribution is computed using eigenvalue decomposition.

## Examples

``` r
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
```
