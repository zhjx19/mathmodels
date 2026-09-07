test_that("HHI works", {
  x = c(50, 30, 20)
  r = HHI(x)
  expect_equal(r, sum((x/sum(x))^2))
  r2 = HHI(x, scaled=TRUE)
  expect_true(r2 >= 0 & r2 <= 1)
})

test_that("HHI single firm returns 1", {
  expect_equal(HHI(10), 1)
  expect_equal(HHI(10, scaled=TRUE), 1)
})

test_that("HHI input validation", {
  expect_error(HHI("a"), "must be a numeric vector")
  expect_error(HHI(c(10, -5, 5)), "negative")
  expect_error(HHI(c(0, 0)), "sum")
})

test_that("EG requires no missing values in y", {
  eg_data = data.frame(
    region = c("R1", "R1", "R2", "R2"),
    industry = c("A", "A", "A", "A"),
    employment = c(250, 200, 150, NA)
  )
  expect_error(EG(eg_data, region, industry, y = employment), "missing values")
})
