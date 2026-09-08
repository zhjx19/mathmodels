# Tests for the three-line table export (export_table.R)

d = data.frame(
  Var = c("mean", "sd", "n"),
  Value = c(20.0906, 6.0269, 32)
)

reg_res = reg_lm(mpg ~ wt + hp, data = mtcars)

test_that("export_table latex output contains booktabs structure", {
  code = export_table(d, output = "latex", caption = "Summary",
    note = "note text")
  expect_type(code, "character")
  expect_true(any(grepl("\\\\toprule", code)))
  expect_true(any(grepl("\\\\bottomrule", code)))
  expect_true(any(grepl("\\\\caption\\{Summary\\}", code)))
  expect_true(any(grepl("footnotesize", code)))
  # underscores escaped
  df2 = data.frame(my_var = 1:2)
  code2 = export_table(df2, output = "latex")
  expect_true(any(grepl("my\\\\_var", code2)))
})

test_that("export_table latex writes to file", {
  tmp = tempfile(fileext = ".tex")
  code = export_table(d, output = "latex", path = tmp)
  expect_true(file.exists(tmp))
  expect_true(length(readLines(tmp)) > 5)
  unlink(tmp)
})

test_that("export_table docx writes a file", {
  tmp = tempfile(fileext = ".docx")
  ft = export_table(d, output = "docx", path = tmp, caption = "Cap")
  expect_true(file.exists(tmp))
  expect_s3_class(ft, "flextable")
  unlink(tmp)
})

test_that("export_table docx requires path", {
  expect_error(export_table(d, output = "docx"), "path is required")
})

test_that("export_table formats regression results", {
  code = export_table(reg_res, output = "latex", caption = "Model 1")
  expect_type(code, "character")
  expect_true(any(grepl("wt", code)))
  expect_true(any(grepl("\\*|\\\\midrule", code)))
  # p values formatted with stars
  expect_true(any(grepl("\\^\\{\\*\\*\\*\\}|\\*", code)))
})

test_that("export_table handles stat_* results", {
  r = stat_t_test(sleep, .cols = extra, group = group)
  code = export_table(r, output = "latex")
  expect_true(any(grepl("\\\\toprule", code)))
  expect_true(any(grepl("p.value", code)))
})

test_that("export_table validates input", {
  expect_error(export_table(list(unrelated = 1), output = "latex"),
    "reg_\\* result")
  expect_error(export_table(3, output = "latex"), NULL)
})
