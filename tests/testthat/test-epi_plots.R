# 传染病可视化与指标函数单元测试 -------------------------------------------

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

# =========================== epi_metrics ====================================

test_that("epi_metrics returns list with 4 summary fields", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  expect_type(r, "list")
  expect_named(r, c("R0", "peak_infection", "peak_time", "attack_rate"))
})

test_that("epi_metrics R0 = beta * N / gamma", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  N = 990 + 10
  expected_R0 = 0.002 * N / 0.1
  expect_equal(r$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epi_metrics uses explicit N", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1, N = 2000)
  expected_R0 = 0.002 * 2000 / 0.1
  expect_equal(r$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epi_metrics peak_infection matches max(I)", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  expect_equal(r$peak_infection, max(sir$I))
})

test_that("epi_metrics attack_rate in [0,1]", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  expect_true(r$attack_rate >= 0 && r$attack_rate <= 1)
})

test_that("epi_metrics infers N from initial state", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  N_inferred = sum(sir[1, setdiff(names(sir), "time")])
  expected_R0 = 0.002 * N_inferred / 0.1
  expect_equal(r$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epi_metrics with SEIR: N inferred from all compartments", {
  r = epi_metrics(seir, beta = 0.3, gamma = 0.1)
  N_inferred = sum(seir[1, setdiff(names(seir), "time")])
  expected_R0 = 0.3 * N_inferred / 0.1
  expect_equal(r$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epi_metrics errors on missing time", {
  df_bad = sir[, setdiff(names(sir), "time"), drop = FALSE]
  expect_error(epi_metrics(df_bad, beta = 0.002, gamma = 0.1), "time")
})

test_that("epi_metrics errors on missing S", {
  df_bad = sir[, c("time", "I"), drop = FALSE]
  expect_error(epi_metrics(df_bad, beta = 0.002, gamma = 0.1), "'S'")
})

test_that("epi_metrics errors on missing I", {
  df_bad = sir[, c("time", "S"), drop = FALSE]
  expect_error(epi_metrics(df_bad, beta = 0.002, gamma = 0.1), "'I'")
})

test_that("epi_metrics errors on missing beta", {
  expect_error(epi_metrics(sir, gamma = 0.1), "Must provide")
})

test_that("epi_metrics errors on missing gamma", {
  expect_error(epi_metrics(sir, beta = 0.002), "Must provide")
})



test_that("plot_incidence returns ggplot object", {
  p = plot_incidence(sir)
  expect_s3_class(p, "ggplot")
})

test_that("plot_incidence works with SEIR data", {
  p = plot_incidence(seir)
  expect_s3_class(p, "ggplot")
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
