test_that("gini0 works", {
  income = c(10, 20, 30, 40, 100)
  r = gini0(income)
  expect_true(r >= 0 & r <= 1)
  expect_equal(gini0(c(10,10,10,10)), 0, tolerance=1e-8)
})

test_that("gini works", {
  y = c(10, 20, 30, 40, 100)
  pop = c(100, 150, 200, 250, 300)
  r = gini(y, pop)
  expect_true(r >= 0 & r <= 1)
})

test_that("theil0 works", {
  r = theil0(c(10, 10, 10))
  expect_equal(r, 0, tolerance=1e-8)
})

test_that("theil works", {
  y = c(10, 20, 30)
  pop = c(5, 5, 5)
  r = theil(y, pop)
  expect_true(r >= 0)
})

test_that("theil0_g works", {
  data = data.frame(g=c("A","A","B","B","B"), y=c(10,10,8,6,6))
  r = theil0_g(data, "g", "y")
  expect_type(r, "list")
  expect_named(r, c("theil", "ratio"))
})

test_that("theil_g works", {
  data = data.frame(g=c("A","A","B"), y=c(10,20,30), pop=c(3,2,5))
  r = theil_g(data, "g", "y", "pop")
  expect_type(r, "list")
})

test_that("inequality input validation", {
  expect_error(gini0("a"), "must be a numeric vector")
  expect_error(gini0(c(-1, 2)), "non-negative")
  expect_error(gini(c(1,2), c(3,4,5)), "same length")
  expect_error(theil_g(1:5, "g", "y", "p"), "must be a data frame")
})
