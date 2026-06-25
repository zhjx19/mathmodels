# 传染病决策指标函数单元测试 -------------------------------------------------

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

# =========================== epidemic_metrics ================================

test_that("epidemic_metrics returns list with summary and trajectory", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_type(r, "list")
  expect_named(r, c("summary", "trajectory"))
})

test_that("epidemic_metrics summary contains all expected fields", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_named(r$summary,
    c("R0", "peak_infection", "peak_time", "attack_rate",
      "final_susceptible", "final_infected", "control_time_steps",
      "time_above_threshold"))
})

test_that("epidemic_metrics trajectory has derived columns", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_true(all(c("Rt", "incidence", "growth_rate") %in% names(r$trajectory)))
})

# =========================== R0 computation ==================================

test_that("epidemic_metrics computes R0 = beta * N / gamma", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  N = 990 + 10  # S+I = 1000
  expected_R0 = 0.002 * N / 0.1  # 0.02
  expect_equal(r$summary$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epidemic_metrics uses explicit N when provided", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1), N = 2000)
  expected_R0 = 0.002 * 2000 / 0.1  # 0.04
  expect_equal(r$summary$R0, expected_R0, ignore_attr = TRUE)
})

# =========================== peak detection ==================================

test_that("epidemic_metrics peak_infection matches max(I)", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_equal(r$summary$peak_infection, max(sir$I))
})

test_that("epidemic_metrics attack_rate in [0,1]", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_true(r$summary$attack_rate >= 0 && r$summary$attack_rate <= 1)
})

# =========================== threshold =======================================

test_that("epidemic_metrics time_above_threshold is NA when threshold=NULL", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_true(is.na(r$summary$time_above_threshold))
})

test_that("epidemic_metrics computes time_above_threshold correctly", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1), threshold = 100)
  manual_count = sum(sir$I > 100, na.rm = TRUE)
  expect_equal(r$summary$time_above_threshold, manual_count)
})

# =========================== control_time ====================================

test_that("epidemic_metrics control_time_steps matches Rt < 1 count", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  manual = sum(r$trajectory$Rt < 1, na.rm = TRUE)
  expect_equal(r$summary$control_time_steps, manual)
})

# =========================== Rt in trajectory ================================

test_that("epidemic_metrics Rt = R0 * S / N at each time", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  N = 990 + 10
  expected_Rt = r$summary$R0 * sir$S / N
  expect_equal(r$trajectory$Rt, expected_Rt, tolerance = 1e-10)
})

# =========================== error handling ==================================

test_that("epidemic_metrics errors on missing time column", {
  df_bad = sir[, setdiff(names(sir), "time"), drop = FALSE]
  expect_error(
    epidemic_metrics(df_bad, params = c(beta = 0.002, gamma = 0.1)),
    "time"
  )
})

test_that("epidemic_metrics errors on missing S column", {
  df_bad = sir[, c("time", "I"), drop = FALSE]
  expect_error(
    epidemic_metrics(df_bad, params = c(beta = 0.002, gamma = 0.1)),
    "'S'"
  )
})

test_that("epidemic_metrics errors on missing I column", {
  df_bad = sir[, c("time", "S"), drop = FALSE]
  expect_error(
    epidemic_metrics(df_bad, params = c(beta = 0.002, gamma = 0.1)),
    "'I'"
  )
})

test_that("epidemic_metrics errors on missing beta", {
  expect_error(
    epidemic_metrics(sir, params = c(gamma = 0.1)),
    "must contain"
  )
})

test_that("epidemic_metrics errors on missing gamma", {
  expect_error(
    epidemic_metrics(sir, params = c(beta = 0.002)),
    "must contain"
  )
})

# =========================== N auto-inference ================================

test_that("epidemic_metrics infers N from initial state", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  N_inferred = sum(sir[1, setdiff(names(sir), "time")])
  expected_R0 = 0.002 * N_inferred / 0.1
  expect_equal(r$summary$R0, expected_R0, ignore_attr = TRUE)
})

test_that("epidemic_metrics with SEIR: N inferred from all compartments", {
  r = epidemic_metrics(seir, params = c(beta = 0.3, gamma = 0.1))
  N_inferred = sum(seir[1, setdiff(names(seir), "time")])
  expected_R0 = 0.3 * N_inferred / 0.1
  expect_equal(r$summary$R0, expected_R0, ignore_attr = TRUE)
})

# =========================== growth_rate =====================================

test_that("epidemic_metrics growth_rate first element is NA", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  expect_true(is.na(r$trajectory$growth_rate[1]))
})

test_that("epidemic_metrics growth_rate computed from diff(log(I))", {
  r = epidemic_metrics(sir, params = c(beta = 0.002, gamma = 0.1))
  # Compare with manual computation
  I_safe = pmax(sir$I, 1e-8)
  manual_growth = c(NA, diff(log(I_safe)))
  expect_equal(r$trajectory$growth_rate, manual_growth)
})
