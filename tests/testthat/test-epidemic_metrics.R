# 传染病指标函数单元测试 -------------------------------------------------

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

# =========================== epi_metrics =====================================

test_that("epi_metrics returns list with 4 fields", {
  r = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  expect_type(r, "list")
  expect_named(r, c("R0", "peak_infection", "peak_time", "attack_rate"))
})

test_that("epi_metrics accepts named data argument", {
  r = epi_metrics(data = sir, beta = 0.002, gamma = 0.1)
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

# =========================== N auto-inference ================================

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

# =========================== error handling ==================================

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
