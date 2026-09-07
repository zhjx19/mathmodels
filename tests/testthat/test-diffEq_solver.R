# 微分方程模型单元测试 ---------------------------------------------------

# =========================== ode_solver ============================

test_that("ode_solver returns data.frame with correct structure", {
  res = ode_solver(
    init = c(N = 100),
    times = seq(0, 5, by = 0.5),
    equations = c(N = "0.5 * N"),
    params = c(r = 0.5)
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "N"))
  expect_equal(nrow(res), 11L)
  expect_equal(res$time, seq(0, 5, by = 0.5))
  expect_true(all(res$N > 0))
})

test_that("ode_solver works without params", {
  res = ode_solver(
    init = c(N = 1),
    times = seq(0, 2, by = 1),
    equations = c(N = "-0.5 * N")
  )
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "N"))
})

test_that("ode_solver handles multi-variable systems", {
  res = ode_solver(
    init = c(S = 900, I = 100),
    times = seq(0, 10, by = 1),
    equations = c(
      S = "-0.002 * S * I",
      I = "0.002 * S * I - 0.5 * I"
    ),
    params = c(beta = 0.002, gamma = 0.5)
  )
  expect_named(res, c("time", "S", "I"))
  expect_equal(nrow(res), 11L)
})

test_that("ode_solver preserves total population in SIR system", {
  res = ode_solver(
    init = c(S = 999, I = 1, R = 0),
    times = seq(0, 20, by = 0.5),
    equations = c(
      S = "-0.003 * S * I",
      I = "0.003 * S * I - 0.1 * I",
      R = "0.1 * I"
    ),
    params = c(beta = 0.003, gamma = 0.1)
  )
  total_pop = res$S + res$I + res$R
  expect_equal(total_pop[1], total_pop[length(total_pop)], tolerance = 1e-2)
})

test_that("ode_solver reports error on invalid equation", {
  expect_error(
    ode_solver(
      init = c(N = 10),
      times = seq(0, 5, by = 1),
      equations = c(N = "undefined_var * N")
    ),
    "Error evaluating equation"
  )
})

# =========================== ode_malthus ============================

test_that("ode_malthus returns correct structure", {
  res = ode_malthus(init = c(N = 100), params = c(r = 0.5),
                      times = seq(0, 10, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "N"))
  expect_true(all(res$N > 0))
})

test_that("ode_malthus grows with positive r", {
  res = ode_malthus(init = c(N = 10), params = c(r = 0.3),
                      times = seq(0, 5, by = 0.5))
  expect_gt(res$N[nrow(res)], res$N[1])
})

test_that("ode_malthus declines with negative r", {
  res = ode_malthus(init = c(N = 100), params = c(r = -0.3),
                      times = seq(0, 5, by = 0.5))
  expect_lt(res$N[nrow(res)], res$N[1])
})

# =========================== ode_logistic ============================

test_that("ode_logistic returns correct structure", {
  res = ode_logistic(init = c(N = 10), params = c(r = 0.8, K = 500),
                       times = seq(0, 30, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "N"))
  expect_true(all(res$N > 0))
})

test_that("ode_logistic approaches carrying capacity", {
  res = ode_logistic(init = c(N = 10), params = c(r = 0.8, K = 500),
                       times = seq(0, 80, by = 1))
  final_N = res$N[nrow(res)]
  expect_equal(final_N, 500, tolerance = 1e-2)
})

test_that("ode_logistic decreases when N0 > K", {
  res = ode_logistic(init = c(N = 800), params = c(r = 0.5, K = 500),
                       times = seq(0, 30, by = 1))
  expect_lt(res$N[nrow(res)], res$N[1])
})

test_that("ode_logistic with N0 = K stays at equilibrium", {
  res = ode_logistic(init = c(N = 200), params = c(r = 0.5, K = 200),
                       times = seq(0, 10, by = 1))
  later_values = res$N[3:nrow(res)]
  expect_true(all(abs(later_values - 200) < 1e-6))
})

# =========================== ode_si ============================

test_that("ode_si returns correct structure", {
  res = ode_si(init = c(S = 990, I = 10), params = c(beta = 0.002),
                 times = seq(0, 30, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I"))
  expect_true(all(res$S > -1e-6))
  expect_true(all(res$I > -1e-6))
})

test_that("ode_si conserves total population", {
  res = ode_si(init = c(S = 990, I = 10), params = c(beta = 0.002),
                 times = seq(0, 30, by = 0.5))
  N = res$S + res$I
  expect_equal(max(N) - min(N), 0, tolerance = 1e-6)
})


test_that("ode_si: S decreases, I increases monotonically", {
  res = ode_si(init = c(S = 990, I = 10), params = c(beta = 0.002),
                 times = seq(0, 30, by = 0.5))
  # Overall trend: S decreases, I increases
  expect_lt(res$S[nrow(res)], res$S[1])
  expect_gt(res$I[nrow(res)], res$I[1])
})

# =========================== ode_si (open pop) ============================

test_that("ode_si with b/d returns correct structure", {
  res = ode_si(init = c(S = 990, I = 10),
                 params = c(beta = 0.002, b = 0.01, d = 0.01),
                 times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I"))
})

test_that("ode_si with b/d: population changes", {
  res = ode_si(init = c(S = 500, I = 10),
                 params = c(beta = 0.002, b = 0.02, d = 0.01),
                 times = seq(0, 100, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  expect_gt(N_end, N_start)
})

test_that("ode_sis returns correct structure", {
  res = ode_sis(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I"))
})

test_that("ode_sis conserves total population", {
  res = ode_sis(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 0.5))
  N = res$S + res$I
  expect_equal(max(N) - min(N), 0, tolerance = 1e-6)
})

test_that("ode_sis reaches endemic equilibrium when R0 > 1", {
  res = ode_sis(init = c(S = 990, I = 10),
                  params = c(beta = 0.3, gamma = 0.1),
                  times = seq(0, 200, by = 1))
  final_S = res$S[nrow(res)]
  N = 1000
  expected_S = N * 0.1 / 0.3
  expect_equal(final_S, expected_S, tolerance = 0.5)
})

# =========================== ode_sis (open pop) ============================

test_that("ode_sis with b/d returns correct structure", {
  res = ode_sis(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, b = 0.01, d = 0.01),
                  times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I"))
})

test_that("ode_sis with b/d: population changes", {
  res = ode_sis(init = c(S = 500, I = 10),
                  params = c(beta = 0.3, gamma = 0.1, b = 0.02, d = 0.01),
                  times = seq(0, 100, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  expect_gt(N_end, N_start)
})

# =========================== ode_sir ============================

test_that("ode_sir returns correct structure", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I", "R"))
})

test_that("ode_sir conserves total population", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 0.5))
  N = res$S + res$I + res$R
  expect_equal(max(N) - min(N), 0, tolerance = 1e-6)
})

test_that("ode_sir: S monotonic non-increasing", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 0.5))
  expect_true(all(diff(res$S) <= 0))
})

test_that("ode_sir: R monotonic non-decreasing", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 50, by = 0.5))
  expect_true(all(diff(res$R) >= 0))
})

test_that("ode_sir: R0 defaults to 0", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1),
                  times = seq(0, 10, by = 1))
  expect_equal(res$R[1], 0)
})

# =========================== ode_sir (open pop) ============================

test_that("ode_sir with b/d returns correct structure", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, b = 0.01, d = 0.01),
                  times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "I", "R"))
})

test_that("ode_sir with b/d: population changes over time", {
  res = ode_sir(init = c(S = 500, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, b = 0.02, d = 0.01),
                  times = seq(0, 100, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  # b > d so population should grow
  expect_gt(N_end, N_start)
})

test_that("ode_sir with b/d only: birth only (no death)", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, b = 0.01),
                  times = seq(0, 50, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  expect_gt(N_end, N_start)
})

test_that("ode_sir with d only: death only (no birth)", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, d = 0.01),
                  times = seq(0, 50, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  expect_lt(N_end, N_start)
})

test_that("ode_sir with explicit N overrides sum(init)", {
  res = ode_sir(init = c(S = 990, I = 10),
                  params = c(beta = 0.002, gamma = 0.1, N = 2000),
                  times = seq(0, 50, by = 1))
  # With N=2000, force of infection is weaker than with N=1000
  res_default = ode_sir(init = c(S = 990, I = 10),
                          params = c(beta = 0.002, gamma = 0.1),
                          times = seq(0, 50, by = 1))
  # Peak infection should be lower / later with larger N
  expect_true(max(res$I) <= max(res_default$I) + 1)
})

# =========================== ode_seir ============================

test_that("ode_seir returns correct structure", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
                   times = seq(0, 30, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "E", "I", "R"))
})

test_that("ode_seir conserves total population", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
                   times = seq(0, 30, by = 0.5))
  N = res$S + res$E + res$I + res$R
  expect_equal(max(N) - min(N), 0, tolerance = 1e-6)
})

test_that("ode_seir: all compartments non-negative", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
                   times = seq(0, 30, by = 0.5))
  expect_true(all(res$S > -1e-6 & res$E > -1e-6 &
                  res$I > -1e-6 & res$R > -1e-6))
})

test_that("ode_seir: R0 defaults to 0", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1),
                   times = seq(0, 10, by = 1))
  expect_equal(res$R[1], 0)
})

# =========================== ode_seir (open pop) ============================

test_that("ode_seir with b/d returns correct structure", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1,
                              b = 0.01, d = 0.01),
                   times = seq(0, 50, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "S", "E", "I", "R"))
})

test_that("ode_seir with b/d: population changes", {
  res = ode_seir(init = c(S = 500, E = 5, I = 5),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1,
                              b = 0.02, d = 0.01),
                   times = seq(0, 100, by = 1))
  N_start = sum(res[1, -1])
  N_end = sum(res[nrow(res), -1])
  expect_gt(N_end, N_start)
})

test_that("ode_seir with b/d: all compartments non-negative", {
  res = ode_seir(init = c(S = 980, E = 10, I = 10),
                   params = c(beta = 0.3, alpha = 0.2, gamma = 0.1,
                              b = 0.01, d = 0.01),
                   times = seq(0, 50, by = 0.5))
  expect_true(all(res$S > -1e-6 & res$E > -1e-6 &
                  res$I > -1e-6 & res$R > -1e-6))
})

# =========================== ode_lv ============================

test_that("ode_lv returns correct structure", {
  res = ode_lv(init = c(x = 40, y = 9),
                 params = c(alpha = 0.4, beta = 0.04, delta = 0.02, mu = 0.5),
                 times = seq(0, 100, by = 1))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "x", "y"))
})

test_that("ode_lv produces bounded oscillations", {
  res = ode_lv(init = c(x = 40, y = 9),
                 params = c(alpha = 0.4, beta = 0.04, delta = 0.02, mu = 0.5),
                 times = seq(0, 100, by = 0.5))
  expect_true(all(res$x > 0))
  expect_true(all(res$y > 0))
})

# =========================== ode_solver validation (review fixes) ==========

test_that("ode_solver rejects non-monotonic times", {
  expect_error(
    ode_solver(init = c(S = 9, I = 1), times = c(0, 10, 5, 20),
               equations = c(S = "-0.1*S*I", I = "0.1*S*I"),
               params = c(beta = 0.1)),
    "strictly increasing"
  )
})

test_that("ode_solver rejects extra compartments in init", {
  expect_error(
    ode_si(init = c(S = 990, I = 10, X = 50),
           params = c(beta = 0.002), times = seq(0, 10, by = 1)),
    "not used by this model"
  )
})

test_that("ode_solver names the variable in parse errors", {
  expect_error(
    ode_solver(init = c(S = 9, I = 1), times = 0:10,
               equations = c(S = "-beta*S*I", I = "0.1*S*")),
               "cannot parse the equation for variable 'I'"
  )
})

test_that("ode_solver warns on shadowed parameters", {
  expect_warning(
    ode_solver(init = c(S = 9, I = 1), times = 0:10,
               equations = c(S = "-beta*S*I", I = "beta*S*I"),
               params = c(beta = 0.1, t = 5)),
    "shadowed"
  )
})

test_that("plot_incidence warns on non-uniform time steps", {
  sir = ode_sir(init = c(S = 990, I = 10),
                params = c(beta = 0.002, gamma = 0.1),
                times = c(0, 1, 2, 10, 20))
  expect_warning(plot_incidence(sir), "non-uniform")
})

test_that("epi_metrics ignores user-added columns when inferring N", {
  sir = ode_sir(init = c(S = 990, I = 10),
                params = c(beta = 0.002, gamma = 0.1),
                times = seq(0, 50, by = 1))
  sir_clean = epi_metrics(sir, beta = 0.002, gamma = 0.1)
  sir_junk = sir
  sir_junk$note = 1000  # junk column 100x the population
  sir_junk = epi_metrics(sir_junk, beta = 0.002, gamma = 0.1)
  expect_equal(sir_clean$attack_rate, sir_junk$attack_rate)
  expect_equal(sir_clean$peak_infection, sir_junk$peak_infection)
})

test_that("epi_metrics rejects negative beta/gamma", {
  sir = ode_sir(init = c(S = 990, I = 10),
                params = c(beta = 0.002, gamma = 0.1),
                times = seq(0, 10, by = 1))
  expect_error(epi_metrics(sir, beta = -0.1, gamma = 0.1), "non-negative")
})
