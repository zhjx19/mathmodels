test_that("rank_sum_ratio works with integer method", {
  data = data.frame(ID=c("A","B","C"), X1=c(10,20,15), X2=c(5,10,8))
  w = c(0.4, 0.6)
  r = rank_sum_ratio(data, w, method="int")
  expect_type(r, "list")
  expect_named(r, c("resultTable", "reg", "rankTable"))
  expect_s3_class(r$reg, "lm")
})

test_that("rank_sum_ratio works with non-int method", {
  data = data.frame(ID=c("A","B","C"), X1=c(10,20,15), X2=c(5,10,8))
  r = rank_sum_ratio(data, method="non-int")
  expect_s3_class(r$reg, "lm")
})

test_that("rank_sum_ratio input validation", {
  data = data.frame(ID=c("A"), X1=c(10))
  expect_error(rank_sum_ratio(1:5), "must be a data frame")
})
