test_that("topsis works", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  w = c(0.6, 0.4)
  idx = c("+","-")
  r = topsis(A, w, idx)
  expect_true(all(r >= 0 & r <= 1))
  expect_length(r, 3)
})

test_that("topsis defaults", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  r = topsis(A)
  expect_true(all(r >= 0 & r <= 1))
})

test_that("topsis input validation", {
  expect_error(topsis(data.frame(x=c(1,2))), "at least 2 columns")
  expect_error(topsis(data.frame(x1=c(2,5,3), x2=c(8,1,6)), w=c(0.5)), "must have length equal to ncol")
})

test_that("topsis rejects invalid index values", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  expect_error(topsis(A, index=c("+","x")), "index must contain only")
})

test_that("topsis rejects all-zero criterion columns", {
  A = data.frame(X1=c(2,5,3), X2=c(0,0,0))
  expect_error(topsis(A), "only zeros")
})

test_that("topsis names samples when X has no rownames", {
  A = matrix(c(2,5,3, 8,1,6), nrow=3)
  r = topsis(A)
  expect_named(r, c("Sample1","Sample2","Sample3"))
  B = A
  rownames(B) = c("a","b","c")
  expect_named(topsis(B), c("a","b","c"))
})

test_that("topsis rejects NA in w", {
  A = data.frame(X1=c(2,5,3), X2=c(8,1,6))
  expect_error(topsis(A, w=c(0.5, NA)), "no NA")
})

test_that("grey_corr identical comparison sequences return correlation 1", {
  ref = c(1, 2, 3, 4)
  cmp = data.frame(x1 = ref, x2 = ref)
  expect_equal(grey_corr(ref, cmp), c(1, 1))
})

test_that("grey_corr rejects NA or negative w", {
  ref = c(1, 2, 3)
  cmp = data.frame(x1 = c(1, 2, 4))
  expect_error(grey_corr(ref, cmp, w = c(0.5, NA, 0.5)), "without NA")
  expect_error(grey_corr(ref, cmp, w = c(-0.5, 1, 0.5)), "w must")
  expect_error(grey_corr(ref, cmp, w = c(0.5, 0.5)), "length")
})
