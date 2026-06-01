test_that("combine_weights linear", {
  w_subj = c(0.4, 0.3, 0.2, 0.1)
  w_obj = c(0.25, 0.2, 0.3, 0.25)
  r = combine_weights(w_subj, w_obj, "linear", alpha=0.6)
  expect_equal(sum(r), 1, tolerance=1e-8)
  expect_equal(r, 0.6*w_subj + 0.4*w_obj)
})

test_that("combine_weights multiplicative", {
  w_subj = c(0.4, 0.3, 0.2, 0.1)
  w_obj = c(0.25, 0.2, 0.3, 0.25)
  r = combine_weights(w_subj, w_obj, "multiplicative")
  expect_equal(sum(r), 1, tolerance=1e-8)
})

test_that("combine_weights game", {
  r = combine_weights(c(0.4,0.3,0.2,0.1), c(0.25,0.2,0.3,0.25), "game")
  expect_equal(sum(r), 1, tolerance=1e-8)
})

test_that("combine_weights game_linear", {
  r = combine_weights(c(0.4,0.3,0.2,0.1), c(0.25,0.2,0.3,0.25), "game_linear")
  expect_equal(sum(r), 1, tolerance=1e-8)
})

test_that("combine_weights input validation", {
  expect_error(combine_weights(c(-1, 1), c(0.5, 0.3)), "non-negative")
  expect_error(combine_weights(c(0.3, 0.3), c(0.5)), "same length")
})
