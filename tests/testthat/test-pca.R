test_that("pca_weight works", {
  ind = c("+","+","-","-")
  r = pca_weight(iris[1:10, 1:4], ind, nfs = 2)
  expect_type(r, "list")
  expect_equal(sum(r$w), 1, tolerance=1e-8)
  expect_length(r$s, 10)
})

test_that("pca_weight input validation", {
  X = data.frame(x1=c(3,5,2,7))
  expect_error(pca_weight(X), "at least 2 columns")
})
