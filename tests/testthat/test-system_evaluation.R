test_that("coupling_degree works", {
  df = data.frame(
    ID = LETTERS[1:3],
    s1 = c(0.5, 0.7, 0.3),
    s2 = c(0.4, 0.6, 0.8)
  )
  r = coupling_degree(df, id_cols="ID")
  expect_s3_class(r, "tbl_df")
  expect_true(all(r$coupling >= 0 & r$coupling <= 1, na.rm=TRUE))
  expect_true(all(r$coupling_coord >= 0 & r$coupling_coord <= 1, na.rm=TRUE))
})

test_that("coupling_degree adjusted", {
  df = data.frame(s1=c(0.5,0.7,0.3), s2=c(0.4,0.6,0.8), s3=c(0.6,0.5,0.4))
  r = coupling_degree(df, type="adjusted")
  expect_true(all(r$coupling >= 0 & r$coupling <= 1, na.rm=TRUE))
})

test_that("obstacle_degree works", {
  df = data.frame(ID=LETTERS[1:3], x1=c(0.8,0.6,0.4), x2=c(0.5,0.7,0.3))
  r = obstacle_degree(df, id_cols="ID")
  expect_s3_class(r, "tbl_df")
  expect_named(r, c("ID", "x1", "x2"))
})

test_that("obstacle_degree scaled", {
  df = data.frame(x1=c(0.8,0.6,0.4), x2=c(0.5,0.7,0.3))
  r = obstacle_degree(df, scaled=TRUE)
  expect_equal(rowSums(as.data.frame(r)), c(1,1,1), tolerance=1e-8)
})

test_that("system_evaluation input validation", {
  expect_error(coupling_degree("a"), "must be a data frame")
  expect_error(obstacle_degree("a"), "must be a data frame")
})

test_that("coupling_degree input validation", {
  df = data.frame(s1=c(0.5,0.7), s2=c(0.4,0.6))
  expect_error(coupling_degree(df, type="wrong"), "should be one of")
  expect_error(coupling_degree(df, w=c(0.5, 0.3)), "sum to 1")
  expect_error(coupling_degree(df, id_cols="IDX"), "id_cols")
  expect_error(obstacle_degree(df, id_cols="IDX"), "id_cols")
})

test_that("obstacle_degree scaled all-one row is neutral 0.5", {
  df = data.frame(x1=c(1, 0.6), x2=c(1, 0.4))
  expect_warning(r <- obstacle_degree(df, scaled=TRUE), "0.5")
  expect_equal(as.numeric(r[1, ]), c(0.5, 0.5))
  expect_equal(sum(as.numeric(r[2, ])), 1, tolerance=1e-8)
})

