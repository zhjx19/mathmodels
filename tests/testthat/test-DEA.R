# ============================================================================
# DEA Models Tests
# ============================================================================

library(testthat)

test_that("basic_DEA returns correct structure", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")

  expect_type(r, "list")
  expect_named(r, c("efficiencies", "lambdas", "slacks", "targets",
                    "returns", "model", "orientation", "dmu"))
  expect_length(r$efficiencies, 7L)
  expect_true(all(r$efficiencies >= 0 & r$efficiencies <= 1 + 1e-10, na.rm = TRUE))
})

test_that("basic_DEA CRS vs VRS", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r_crs = basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
  r_vrs = basic_DEA(df, inputs = 2:3, outputs = 4, rts = "vrs")

  # VRS efficiencies should be >= CRS efficiencies (larger feasible region)
  expect_true(all(r_vrs$efficiencies >= r_crs$efficiencies - 1e-10))
})

test_that("basic_DEA output orientation returns consistent values", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = basic_DEA(df, inputs = 2:3, outputs = 4, orientation = "oo", rts = "crs")

  expect_type(r, "list")
  expect_length(r$efficiencies, 7L)
  expect_true(all(r$efficiencies >= 0 & r$efficiencies <= 1 + 1e-10, na.rm = TRUE))
})

test_that("super_DEA excludes self", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = super_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")

  expect_length(r$efficiencies, 7L)
  # Super-efficiency can be > 1
  expect_true(any(r$efficiencies > 1, na.rm = TRUE))
  # DMU7 is the most efficient (lowest inputs, moderate output)
  expect_true(!is.na(r$efficiencies[7]))
})

test_that("basic_SBM returns valid efficiency values", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = basic_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")

  expect_length(r$efficiencies, 7L)
  expect_true(all(r$efficiencies >= 0 & r$efficiencies <= 1 + 1e-10, na.rm = TRUE))
  # SBM efficiencies should be <= radial efficiencies (SBM is stricter)
  r_ccr = basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
  expect_true(all(r$efficiencies <= r_ccr$efficiencies + 1e-10))
})

test_that("super_SBM handles super-efficient DMUs gracefully", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = super_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")

  expect_length(r$efficiencies, 7L)
  # DMU1 and DMU7 may be NA if LP infeasible
  expect_true(is.na(r$efficiencies[1]) || r$efficiencies[1] >= 1)
})

test_that("malmquist returns correct structure", {
  panel = data.frame(
    DMU    = rep(paste0("DMU", 1:5), 3),
    Period = rep(1:3, each = 5),
    x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
    y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
  )

  r = malmquist(panel, period = "Period", inputs = 3, outputs = 4,
    rts = "crs", type1 = "cont", type2 = "fgnz")

  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 10L)  # (3-1) * 5
  expect_named(r, c("Period", "DMU", "mi", "ec", "tc", "pech", "sech"))
  expect_true(all(r$mi > 0, na.rm = TRUE))
})

test_that("malmquist supports all type1 × type2 combos", {
  panel = data.frame(
    DMU    = rep(paste0("DMU", 1:5), 3),
    Period = rep(1:3, each = 5),
    x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
    y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
  )

  for (t1 in c("cont", "seq", "glob")) {
    for (t2 in c("fgnz", "rd")) {
      r = malmquist(panel, period = "Period", inputs = 3, outputs = 4,
        rts = "crs", type1 = t1, type2 = t2)
      expect_equal(nrow(r), 10L)
    }
  }
})

test_that("DEA with undesirable outputs works", {
  df = data.frame(
    DMU  = paste0("DMU", 1:5),
    x1   = c(4, 7, 8, 4, 2),
    x2   = c(3, 5, 2, 8, 7),
    y_g  = c(2, 5, 6, 3, 1),
    y_b  = c(1, 2, 2.5, 1.5, 0.5)
  )

  r = basic_DEA(df, inputs = 2:3, outputs = 4:5, ud_outputs = 2, rts = "crs")

  expect_length(r$efficiencies, 5L)
  # Slack matrix for undesirable outputs should exist
  expect_equal(ncol(r$slacks$outputs), 2L)
})

test_that("DEA validates inputs correctly", {
  df = data.frame(
    DMU = paste0("DMU", 1:3),
    x1  = c(1, NA, 3),
    y1  = c(1, 2, 3)
  )
  expect_error(basic_DEA(df, inputs = 2, outputs = 3), "NA")

  # Fewer than 3 columns
  expect_error(basic_DEA(df[, 1:2], inputs = 2, outputs = 2), "at least 3")
})
