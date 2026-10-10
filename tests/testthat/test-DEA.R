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

test_that("super_SBM implements Tone (2002) super-efficiency", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )

  r = super_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")

  expect_length(r$efficiencies, 7L)
  # No NA: efficient DMUs get rho >= 1, inefficient ones get exactly 1
  expect_true(all(!is.na(r$efficiencies)))
  expect_true(all(r$efficiencies >= 1 - 1e-10))
  # Inefficient DMUs (basic_SBM < 1) must have super-SBM exactly 1
  sbm = basic_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
  ineff = sbm$efficiencies < 1 - 1e-10
  expect_true(all(abs(r$efficiencies[ineff] - 1) < 1e-8))
})

test_that("super_SBM hand-check: Tone (2002) example", {
  # 1 input, 1 output, CRS: A(1,1) efficient, B(2,0.5).
  # Super-efficiency of A against {B}: rho = 4 (hand-derived in review).
  df = data.frame(
    DMU = c("A", "B"),
    x   = c(1, 2),
    y   = c(1, 0.5)
  )
  r = super_SBM(df, inputs = 2, outputs = 3, rts = "crs")
  expect_equal(unname(r$efficiencies["A"]), 4, tolerance = 1e-8)
  # B is inefficient: super-SBM = 1
  expect_equal(unname(r$efficiencies["B"]), 1, tolerance = 1e-8)
})

test_that("SBM with undesirable outputs uses Tone (2004) weights", {
  # CRS, A(x=1, yg=1, yb=1), B(x=2, yg=1, yb=2): SBM-UD score of B = 1/3.
  # Derivation: rho = (1 - s-/x0) / (1 + sg/yg0 + sb/(2*yb0/2)) with
  # x0 = 2 - s- = lambda, yg0 = 1 = lambda - sg, yb0 = 2 = lambda + sb;
  # optimum at lambda = 1: s- = 1, sg = 0, sb = 1 ->
  # rho = (1 - 1/2) / (1 + 0/1 + 1/2) = (1/2)/(3/2) = 1/3.
  df = data.frame(
    DMU = c("A", "B"),
    x   = c(1, 2),
    yg  = c(1, 1),
    yb  = c(1, 2)
  )
  r = basic_SBM(df, inputs = 2, outputs = 3:4, ud_outputs = 2, rts = "crs")
  expect_equal(unname(r$efficiencies["B"]), 1 / 3, tolerance = 1e-8)
  # Targets now include undesirable outputs (Y_b - s_b)
  expect_equal(nrow(r$targets$outputs), 2L)
})

test_that("basic_DEA ud_outputs changes efficiencies (treat-as-input)", {
  df = data.frame(
    DMU  = paste0("DMU", 1:5),
    x1   = c(4, 7, 8, 4, 2),
    x2   = c(3, 5, 2, 8, 7),
    y_g  = c(2, 5, 6, 3, 1),
    y_b  = c(1, 2, 2.5, 1.5, 0.5)
  )
  r_no = basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
  r_ud = basic_DEA(df, inputs = 2:3, outputs = 4:5, ud_outputs = 2, rts = "crs")
  expect_false(isTRUE(all.equal(r_no$efficiencies, r_ud$efficiencies)))
  # UD slacks are reported (not NA)
  expect_true(all(!is.na(r_ud$slacks$outputs[, 2])))
})

test_that("SBM rejects orientation argument", {
  df = data.frame(DMU = paste0("DMU", 1:3), x = c(1, 2, 3), y = c(1, 1, 2))
  expect_error(basic_SBM(df, inputs = 2, outputs = 3, orientation = "oo"),
               "non-oriented")
  expect_error(super_SBM(df, inputs = 2, outputs = 3, orientation = "oo"),
               "non-oriented")
})

test_that("DEA rejects non-positive data", {
  df = data.frame(DMU = paste0("DMU", 1:3), x = c(1, 0, 3), y = c(1, 1, 2))
  expect_error(basic_DEA(df, inputs = 2, outputs = 3), "strictly positive")
  expect_error(basic_SBM(df, inputs = 2, outputs = 3), "strictly positive")
  df2 = data.frame(DMU = paste0("DMU", 1:3), x = c(1, 2, 3), y = c(1, 1, 2))
  panel_bad = rbind(
    data.frame(DMU = df2$DMU, Period = 1, x = df2$x, y = df2$y),
    data.frame(DMU = df2$DMU, Period = 2, x = c(1, 0, 3), y = df2$y)
  )
  expect_error(malmquist(panel_bad, period = "Period", inputs = 3, outputs = 4,
                         rts = "crs"), "strictly positive")
})

test_that("malmquist sequential differs from contemporaneous", {
  panel = data.frame(
    DMU    = rep(paste0("DMU", 1:5), 3),
    Period = rep(1:3, each = 5),
    x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
    y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
  )
  r_cont = malmquist(panel, period = "Period", inputs = 3, outputs = 4,
    rts = "crs", type1 = "cont", type2 = "fgnz")
  r_seq = malmquist(panel, period = "Period", inputs = 3, outputs = 4,
    rts = "crs", type1 = "seq", type2 = "fgnz")
  # nt = 3: transitions 2->3 use the cumulative frontier, so seq differs
  expect_false(isTRUE(all.equal(r_cont$tc, r_seq$tc)))
  # mi consistency: mi = ec * tc for seq
  expect_equal(r_seq$mi, r_seq$ec * r_seq$tc, tolerance = 1e-10)
})

test_that("malmquist sorts numeric-looking character periods numerically", {
  panel_chr = data.frame(
    DMU    = rep(paste0("DMU", 1:3), 2),
    Period = rep(c("2", "10"), each = 3),
    x1     = c(10, 20, 15, 12, 22, 17),
    y1     = c(100, 150, 120, 110, 160, 130)
  )
  panel_num = panel_chr
  panel_num$Period = rep(c(2, 10), each = 3)

  r_chr = malmquist(panel_chr, period = "Period", inputs = 3, outputs = 4,
    rts = "crs", type1 = "cont", type2 = "fgnz")
  r_num = malmquist(panel_num, period = "Period", inputs = 3, outputs = 4,
    rts = "crs", type1 = "cont", type2 = "fgnz")
  # Numeric sorting: character periods must give identical results to numeric
  expect_equal(as.data.frame(r_chr), as.data.frame(r_num))
  expect_true(all(r_chr$Period == "2~10"))  # not "10~2"
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

test_that("basic_DEA and super_DEA stay silent on a well-posed problem", {
  df = data.frame(
    DMU = paste0("DMU", 1:7),
    x1  = c(20, 60, 40, 60, 70, 30, 50),
    x2  = c(151, 200, 120, 170, 250, 210, 90),
    y1  = c(100, 210, 150, 240, 220, 80, 200)
  )
  expect_no_warning(basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs"))
  expect_no_warning(super_DEA(df, inputs = 2:3, outputs = 4, rts = "crs"))
})

test_that("super_DEA reports LP-infeasible DMUs (VRS + super-efficient DMU)", {
  # F is super-efficient. Under VRS, once F is removed from the reference set
  # no convex combination of the remaining DMUs can reach its output levels,
  # so its LP is infeasible -- the textbook VRS super-efficiency failure. It
  # needs VRS *and* an extreme DMU: CRS never fails here, which is why random
  # or single-output data does not surface it.
  dea = data.frame(dmu = LETTERS[1:6],
                   x1 = c(2, 3, 4, 5, 3, 4), x2 = c(3, 4, 5, 6, 4, 5),
                   y1 = c(8, 9, 10, 12, 11, 13), y2 = c(5, 6, 7, 8, 7, 9))

  # input-oriented: F fails
  expect_warning(r_io <- super_DEA(dea, 2:3, 4:5, orientation = "io", rts = "vrs"),
                 "LP infeasible for DMU")
  expect_true(is.na(r_io$efficiencies[["F"]]))
  expect_false(anyNA(r_io$efficiencies[setdiff(names(r_io$efficiencies), "F")]))

  # output-oriented: A fails
  expect_warning(r_oo <- super_DEA(dea, 2:3, 4:5, orientation = "oo", rts = "vrs"),
                 "LP infeasible for DMU")
  expect_true(is.na(r_oo$efficiencies[["A"]]))

  # CRS is feasible on the same data
  expect_no_warning(super_DEA(dea, 2:3, 4:5, orientation = "io", rts = "crs"))
  expect_no_warning(super_DEA(dea, 2:3, 4:5, orientation = "oo", rts = "crs"))
})

test_that("basic_DEA's LP-failure branch is wired (unreachable on real data)", {
  # basic_DEA always keeps the target DMU in its own reference set, so
  # theta = 1 with lambda = e_i is feasible and the solver cannot fail. Unlike
  # super_DEA the branch is therefore unreachable from real data and has to be
  # forced.
  local_mocked_bindings(
    solve_radial_lp = function(...) {
      list(theta = NA_real_, phi = NA_real_, lambda = rep(NA_real_, 1),
           s_in = NA_real_, s_out = NA_real_, s_b = numeric(0))
    },
    .package = "mathmodels"
  )
  df = data.frame(DMU = c("A", "B", "C"), x1 = c(1, 2, 3), y1 = c(3, 2, 1))
  expect_warning(r <- basic_DEA(df, inputs = 2, outputs = 3, rts = "vrs"),
                 "basic_DEA: LP infeasible for DMU")
  expect_true(all(is.na(r$efficiencies)))
})
