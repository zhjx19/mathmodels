# 传染病可视化函数单元测试 -------------------------------------------------

sir = model_sir(
  init   = c(S = 990, I = 10),
  params = c(beta = 0.002, gamma = 0.1),
  times  = seq(0, 50, by = 0.5)
)

seir = model_seir(
  init   = c(S = 980, E = 10, I = 10),
  params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
  times  = seq(0, 30, by = 0.5)
)

# =========================== plot_compartments ============================

test_that("plot_compartments returns ggplot with default vars", {
  p = plot_compartments(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_compartments works with explicit compartments", {
  p = plot_compartments(sir, compartments = c("S", "I"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_compartments works with single compartment", {
  p = plot_compartments(sir, compartments = "I")
  expect_s3_class(p, "ggplot")
})

test_that("plot_compartments with SEIR data", {
  p = plot_compartments(seir, compartments = c("S", "E", "I", "R"))
  expect_s3_class(p, "ggplot")
})

# =========================== compute_incidence ============================

test_that("compute_incidence returns data.frame with expected columns", {
  inc = compute_incidence(sir)
  expect_s3_class(inc, "data.frame")
  expect_named(inc, c("time", "new_infection", "new_infection_I"))
})

test_that("compute_incidence first row is NA", {
  inc = compute_incidence(sir)
  expect_true(is.na(inc$new_infection[1]))
  expect_true(is.na(inc$new_infection_I[1]))
})

test_that("compute_incidence new_infection is non-negative", {
  inc = compute_incidence(sir)
  expect_true(all(inc$new_infection >= 0, na.rm = TRUE))
})

test_that("compute_incidence new_infection matches S depletion", {
  inc = compute_incidence(sir)
  # Total new infections from S depletion
  cum_infections = sum(inc$new_infection, na.rm = TRUE)
  expected_total = sir$S[1] - sir$S[nrow(sir)]
  expect_equal(cum_infections, expected_total, tolerance = 1e-6)
})



test_that("plot_incidence returns ggplot object", {
  p = plot_incidence(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_incidence works with SEIR data", {
  p = plot_incidence(seir)
  expect_s3_class(p, "ggplot")
})

# =========================== plot_infectious_curve ============================

test_that("plot_infectious_curve returns ggplot", {
  p = plot_infectious_curve(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_infectious_curve rejects data without I column", {
  df = data.frame(time = 1:10, S = 10:1)
  expect_error(plot_infectious_curve(df), "Column 'I' not found")
})

# =========================== plot_cumulative_infection ============================

test_that("plot_cumulative_infection returns ggplot", {
  p = plot_cumulative_infection(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cumulative_infection rejects data without S column", {
  df = data.frame(time = 1:10, I = 1:10)
  expect_error(plot_cumulative_infection(df), "Column 'S' not found")
})

# =========================== plot_phase_si ============================

test_that("plot_phase_si returns ggplot", {
  p = plot_phase_si(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_phase_si rejects data without S or I", {
  df = data.frame(time = 1:10, R = 1:10)
  expect_error(plot_phase_si(df), "Requires columns S and I")
})

# =========================== plot_Rt_estimate ============================

test_that("plot_Rt_estimate returns ggplot (mechanistic)", {
  p = plot_Rt_estimate(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_s3_class(p, "ggplot")
})

test_that("plot_Rt_estimate returns ggplot (normalized)", {
  p = plot_Rt_estimate(sir, params = c(beta = 0.002, gamma = 0.1),
                       method = "normalized")
  expect_s3_class(p, "ggplot")
})

test_that("plot_Rt_estimate rejects missing params", {
  expect_error(
    plot_Rt_estimate(sir, params = c(beta = 0.002)),
    "params must contain beta and gamma"
  )
})

test_that("plot_Rt_estimate rejects data without S", {
  df = data.frame(time = 1:10, I = 1:10)
  expect_error(
    plot_Rt_estimate(df, params = c(beta = 0.002, gamma = 0.1)),
    "Column 'S' required"
  )
})

test_that("plot_Rt_estimate works with explicit N", {
  p = plot_Rt_estimate(sir, params = c(beta = 0.002, gamma = 0.1), N = 1000)
  expect_s3_class(p, "ggplot")
})
