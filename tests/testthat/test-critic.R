test_that("weight_critic works", {
  X = data.frame(
    x1 = c(3, 5, 2, 7),
    x2 = c(10, 20, 15, 25)
  )
  index = c("+", "-")
  r = weight_critic(X, index)
  expect_type(r, "list")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 4)
})

test_that("weight_critic entropy method", {
  X = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  index = c("+", "-")
  r = weight_critic(X, index, method="entropy")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
})

test_that("weight_critic input validation", {
  X = data.frame(x1 = c(3, 5, 2, 7))
  expect_error(weight_critic(X), "at least 2 columns")
  expect_error(weight_critic(data.frame(x1=c(1,2), x2=c(3,4))), "at least 3 rows")
  expect_error(weight_critic(data.frame(x1=c(1,2,NA), x2=c(3,4,5))), "NA")
})

test_that("weight_critic matrix input equals data.frame input", {
  Xdf = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  index = c("+", "-")
  r1 = weight_critic(Xdf, index)
  r2 = weight_critic(as.matrix(Xdf), index)
  expect_equal(r1$w, r2$w)
  expect_equal(r1$s, r2$s)
  r3 = weight_critic(Xdf, index, method = "entropy")
  r4 = weight_critic(as.matrix(Xdf), index, method = "entropy")
  expect_equal(r3$w, r4$w)
  expect_equal(r3$s, r4$s)
})

test_that("weight_critic index validation", {
  X = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
  expect_error(weight_critic(X, c("+", "pos")), 'index must contain only')
})

test_that("weight_critic entropy negative data errors", {
  X = data.frame(x1=c(3,-5,2,7), x2=c(10,20,15,25))
  expect_error(weight_critic(X, method="entropy"), "non-negative")
})

