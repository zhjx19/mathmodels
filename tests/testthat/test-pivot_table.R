# Tests for pivot_table

test_that("pivot_table aggregates and spreads", {
  d = data.frame(
    region = c("n", "n", "s", "s"),
    year = c(2023, 2024, 2023, 2024),
    sales = c(10, 20, 30, 40),
    profit = c(1, 2, 3, 4)
  )
  r = pivot_table(d, index = region, columns = year, values = sales)
  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 2)
  expect_equal(names(r), c("region", "sales_2023", "sales_2024"))
  expect_equal(r$sales_2023, c(10, 30))
  expect_equal(r$sales_2024, c(20, 40))
})

test_that("pivot_table multiple values get names_sep suffixes", {
  d = data.frame(
    g = rep(c("a", "b"), each = 4),
    h = rep(c("x", "y"), 4),
    v1 = 1:8,
    v2 = 8:1
  )
  r = pivot_table(d, index = g, columns = h, values = c(v1, v2))
  expect_true(all(c("v1_x", "v1_y", "v2_x", "v2_y") %in% names(r)))
  # custom separator
  r2 = pivot_table(d, index = g, columns = h, values = v1, names_sep = ".")
  expect_true(all(c("v1.x", "v1.y") %in% names(r2)))
})

test_that("pivot_table aggfunc and fill_value work", {
  d = data.frame(
    g = c("a", "a", "b"),
    h = c("x", "y", "x"),
    v = c(1, 3, 10)
  )
  r = pivot_table(d, index = g, columns = h, values = v, aggfunc = sum)
  expect_equal(r$v_y[r$g == "a"], 3)
  expect_true(is.na(r$v_y[r$g == "b"]))
  r2 = pivot_table(d, index = g, columns = h, values = v,
    aggfunc = sum, fill_value = 0)
  expect_equal(r2$v_y[r$g == "b"], 0)
})

test_that("pivot_table multiple index columns", {
  d = data.frame(
    a = rep(c("a1", "a2"), each = 4),
    b = rep(c("b1", "b2"), times = 4),
    c = rep(c("c1", "c2"), each = 2, times = 2),
    v = 1:8
  )
  r = pivot_table(d, index = c(a, b), columns = c, values = v,
    aggfunc = sum)
  expect_equal(nrow(r), 4)
  expect_equal(sum(r$v_c1), 1 + 2 + 5 + 6)
})
