test_that("GM11 works", {
  x = c(100, 120, 145, 175, 210)
  r = GM11(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u", "lambda", "rng"))
  expect_length(r$fitted, length(x))
  expect_type(r$pnext, "double")
  expect_type(r$f, "closure")
})

test_that("GM11 input validation", {
  expect_error(GM11("a"), "must be a numeric vector")
  expect_error(GM11(c(1,2,3)), "length >= 4")
})

test_that("DGM21 works", {
  x = c(2.874, 3.278, 3.39, 3.679, 3.77, 3.8)
  r = DGM21(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u"))
  expect_length(r$fitted, length(x))
})

test_that("verhulst works", {
  x = c(4.93, 2.33, 3.87, 4.35, 6.63, 7.15, 5.37, 6.39, 7.81, 8.35)
  r = verhulst(x)
  expect_type(r, "list")
  expect_named(r, c("fitted", "pnext", "f", "mat", "u"))
  expect_length(r$fitted, length(x))
})

test_that("GM1N works", {
  data = data.frame(
    factor1 = c(50, 55, 60, 65, 70),
    factor2 = c(20, 22, 25, 28, 30),
    output = c(100, 120, 145, 175, 210)
  )
  r = GM1N(data)
  expect_type(r, "list")
  expect_length(r$fitted, nrow(data))
})

test_that("GM1N input validation", {
  expect_error(GM1N("a"), "must be a data frame or matrix")
})
