test_that("critic_weight works", {
  X = data.frame(
    x1 = c(3, 5, 2, 7),
    x2 = c(10, 20, 15, 25)
  )
  index = c("+", "-")
  r = critic_weight(X, index)
  expect_type(r, "list")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 4)
})

test_that("critic_weight entropy method", {
  X = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  index = c("+", "-")
  r = critic_weight(X, index, method="entropy")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
})

test_that("critic_weight input validation", {
  X = data.frame(x1 = c(3, 5, 2, 7))
  expect_error(critic_weight(X), "at least 2 columns")
})
