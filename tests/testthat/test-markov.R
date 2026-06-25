test_that("markov_chain works with factor input", {
  S = factor(c(1, 1, 2, 3, 2, 1, 3, 2, 1, 2,
               3, 1, 2, 3, 1, 2, 1, 3, 3, 1,
               3, 3, 2, 1, 1, 3, 2, 2, 1, 2,
               1, 3, 2, 1, 1, 2, 2, 3, 1, 2),
             levels = 1:3, labels = c("低", "中", "高"))
  
  r = markov_chain(S, s0 = "中", n_steps = 5)
  
  expect_type(r, "list")
  expect_named(r, c("trans_mat", "pred_probs", "pred_states", "pi_final"))
  expect_s3_class(r$trans_mat, "table")
  expect_true(is.matrix(r$pred_probs))
  expect_equal(nrow(r$pred_probs), 5)
  expect_equal(ncol(r$pred_probs), 3)
  expect_length(r$pred_states, 5)
  expect_true(all(r$pred_states %in% c("低", "中", "高")))
  expect_length(r$pi_final, 3)
  expect_equal(sum(r$pi_final), 1)
})

test_that("markov_chain works with character vector", {
  S = c("A", "A", "B", "A", "B", "B", "A", "B")
  r = markov_chain(S, s0 = "A", n_steps = 3)
  
  expect_type(r, "list")
  expect_equal(nrow(r$pred_probs), 3)
  expect_equal(ncol(r$pred_probs), 2)
  expect_length(r$pred_states, 3)
})

test_that("markov_chain input validation", {
  S = factor(c("A", "B", "A", "B", "A"))
  
  expect_error(markov_chain(1:3, "A", 5), "must be one of the state levels")
  expect_error(markov_chain(S, "C", 5), "must be one of the state levels")
  expect_error(markov_chain(S, "A", 0), "n_steps must be a positive integer")
  expect_error(markov_chain(S, "A", 1.5), "n_steps must be a positive integer")
  expect_error(markov_chain(S, "A", -1), "n_steps must be a positive integer")
  expect_error(markov_chain(S, 1, 5), "s0 must be a single character string")
  expect_error(markov_chain(1, "A", 5), "S must have length >= 2")
})

test_that("GM11_markov works", {
  X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
  r = GM11_markov(X, n_ahead = 3)
  
  expect_s3_class(r, "data.frame")
  expect_named(r, c("Period", "Raw", "GM11_fitted", "err_state", 
                    "adj_eff", "Markov_adj"))
  expect_equal(nrow(r), length(X) + 3)
  
  # Historical data should have original values
  expect_equal(r$Raw[1:length(X)], X)
  # Future periods should have NA for Raw
  expect_true(all(is.na(r$Raw[(length(X)+1):nrow(r)])))
})

test_that("GM11_markov with custom breaks", {
  X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
  r = GM11_markov(X, n_ahead = 2, breaks = c(-0.05, 0, 0.05))
  
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), length(X) + 2)
})

test_that("markov_chain edge cases", {
  S = factor(c("A", "B", "A", "B", "A"))

  # Transition matrix rows must sum to 1
  r = markov_chain(S, s0 = "A", n_steps = 5)
  expect_equal(unname(rowSums(r$trans_mat)), rep(1, nrow(r$trans_mat)),
               tolerance = 1e-8)

  # Stationary distribution must sum to 1
  expect_equal(sum(r$pi_final), 1, tolerance = 1e-8)

  # Predicted probabilities each row must sum to 1
  prob_sums = unname(rowSums(r$pred_probs))
  expect_equal(prob_sums, rep(1, nrow(r$pred_probs)), tolerance = 1e-8)

  # Matrix input should error
  expect_error(markov_chain(matrix(1:6, nrow = 2), "A", 3),
               "S must be a vector or factor")
})

test_that("markov_chain with single-step", {
  S = factor(c("A", "B", "A", "B", "A", "B"))
  r = markov_chain(S, s0 = "B", n_steps = 1)

  expect_equal(nrow(r$pred_probs), 1L)
  expect_length(r$pred_states, 1L)
})

test_that("GM11_markov additional edge cases", {
  X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)

  # n_ahead = 1
  r = GM11_markov(X, n_ahead = 1)
  expect_equal(length(r$err_state), length(X) + 1)

  # All NA in future Raw
  future_rows = (length(X) + 1):nrow(r)
  expect_true(all(is.na(r$Raw[future_rows])))

  # Historical Raw must equal input X
  expect_equal(r$Raw[1:length(X)], X)
})

test_that("GM11_markov input validation", {
  X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
  
  expect_error(GM11_markov("a"), "X must be a numeric vector")
  expect_error(GM11_markov(c(1, 2, 3)), "requires at least 4 data points")
  expect_error(GM11_markov(X, n_ahead = 0), "n_ahead must be a positive integer")
  expect_error(GM11_markov(X, n_ahead = 1.5), "n_ahead must be a positive integer")
  expect_error(GM11_markov(X, breaks = "a"), "breaks must be a numeric vector")
})
