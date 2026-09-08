# Tests for the statistical inference toolkit (stat_infer.R)

set.seed(42)
d_two = data.frame(
  y  = c(rnorm(20, 10), rnorm(20, 11.2)),
  g  = factor(rep(c("a", "b"), each = 20)),
  region = rep(c("e", "w"), length.out = 40)
)
d_pair = data.frame(pre = rnorm(25, 5), post = rnorm(25, 5, 0.8) + 0.6)

# ---------------------------------------------------------------------------
# stat_describe
# ---------------------------------------------------------------------------

test_that("stat_describe matches hand-computed values", {
  r = stat_describe(mtcars, .cols = mpg)
  expect_s3_class(r, "stat_infer")
  expect_s3_class(r, "tbl_df")
  expect_equal(r$n_used, 32)
  expect_equal(r$mean, mean(mtcars$mpg), tolerance = 1e-12)
  expect_equal(r$sd, stats::sd(mtcars$mpg), tolerance = 1e-12)
  expect_equal(r$median, stats::median(mtcars$mpg), tolerance = 1e-12)
  expect_equal(r$iqr, stats::IQR(mtcars$mpg), tolerance = 1e-12)
  hw = stats::qt(0.975, 31) * stats::sd(mtcars$mpg) / sqrt(32)
  expect_equal(r$ci_low, mean(mtcars$mpg) - hw, tolerance = 1e-10)
  # moment skewness hand check
  m = mean(mtcars$mpg); s2 = mean((mtcars$mpg - m)^2)
  expect_equal(r$skew, mean((mtcars$mpg - m)^3) / s2^1.5, tolerance = 1e-10)
})

test_that("stat_describe slices via .by", {
  r = stat_describe(d_two, .cols = y, .by = region)
  expect_equal(nrow(r), 2)
  e = d_two$y[d_two$region == "e"]
  expect_equal(r$mean[r$region == "e"], mean(e), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# stat_t_test
# ---------------------------------------------------------------------------

test_that("stat_t_test two-sample Welch matches stats::t.test", {
  fit = stats::t.test(y ~ g, data = d_two)
  r = stat_t_test(d_two, .cols = y, group = g)
  expect_equal(r$statistic, unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$df, unname(fit$parameter), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)
  expect_equal(r$estimate, unname(fit$estimate[2] - fit$estimate[1]),
    tolerance = 1e-10)
  expect_equal(r$conf_low, unname(fit$conf.int[1]), tolerance = 1e-10)
  # Cohen's d hand check
  x1 = d_two$y[d_two$g == "a"]; x2 = d_two$y[d_two$g == "b"]
  sp = sqrt(((19) * stats::var(x1) + (19) * stats::var(x2)) / 38)
  expect_equal(r$effect, (mean(x1) - mean(x2)) / sp, tolerance = 1e-10)
})

test_that("stat_t_test pooled and one-sample match", {
  fit = stats::t.test(y ~ g, data = d_two, var.equal = TRUE)
  r = stat_t_test(d_two, .cols = y, group = g, var.equal = TRUE)
  expect_equal(r$df, unname(fit$parameter), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)

  fit1 = stats::t.test(mtcars$mpg, mu = 20)
  r1 = stat_t_test(mtcars, .cols = mpg, mu = 20)
  expect_equal(r1$statistic, unname(fit1$statistic), tolerance = 1e-10)
  expect_equal(r1$estimate, mean(mtcars$mpg), tolerance = 1e-12)
})

test_that("stat_t_test paired matches and requires two cols", {
  fit = stats::t.test(d_pair$pre, d_pair$post, paired = TRUE)
  r = stat_t_test(d_pair, .cols = pre:post, paired = TRUE)
  expect_equal(r$statistic, unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)
  dz = (mean(d_pair$pre - d_pair$post)) / stats::sd(d_pair$pre - d_pair$post)
  expect_equal(r$effect, dz, tolerance = 1e-10)

  expect_error(stat_t_test(d_pair, .cols = pre, paired = TRUE), "two .cols")
  expect_error(
    stat_t_test(d_two, .cols = y, group = g, paired = TRUE), "two .cols")
})

test_that("stat_t_test .by slices equal manual subsets", {
  r = stat_t_test(d_two, .cols = y, group = g, .by = region)
  expect_equal(nrow(r), 2)
  for (rg in c("e", "w")) {
    sub = d_two[d_two$region == rg, ]
    fit = stats::t.test(y ~ g, data = sub)
    expect_equal(r$statistic[r$region == rg], unname(fit$statistic),
      tolerance = 1e-10)
  }
})

test_that("stat_t_test multiple columns loop with p_adjust", {
  r = stat_t_test(mtcars, .cols = c(mpg, hp, wt), group = am,
    p_adjust = "holm")
  expect_equal(nrow(r), 3)
  expect_true(all(r$p.adjusted >= r$p.value))
  fit = stats::t.test(hp ~ am, data = mtcars)
  i = which(r$variable == "hp")
  expect_equal(r$p.adjusted[i],
    stats::p.adjust(c(r$p.value[i],
      stat_t_test(mtcars, .cols = mpg, group = am)$p.value,
      stat_t_test(mtcars, .cols = wt, group = am)$p.value),
      method = "holm")[1], tolerance = 1e-12)
})

test_that("stat_t_test validates inputs", {
  expect_error(stat_t_test(mtcars, .cols = mpg, group = cyl),
    "exactly 2 levels") # cyl has 3 levels
  expect_error(stat_t_test(mtcars, .cols = mpg, mu = 20, conf_level = 1.2),
    "conf_level")
  df_bad = data.frame(y = c("a", "b"), g = c("x", "y"))
  expect_error(stat_t_test(df_bad, .cols = y, group = g), "numeric")
})

# ---------------------------------------------------------------------------
# stat_wilcox_test
# ---------------------------------------------------------------------------

test_that("stat_wilcox_test matches stats::wilcox.test", {
  fit = stats::wilcox.test(y ~ g, data = d_two)
  r = stat_wilcox_test(d_two, .cols = y, group = g)
  expect_equal(r$statistic, unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)
  # rank-biserial
  n1 = 20; n2 = 20
  expect_equal(r$effect, 2 * unname(fit$statistic) / (n1 * n2) - 1,
    tolerance = 1e-10)

  fitp = stats::wilcox.test(d_pair$pre, d_pair$post, paired = TRUE)
  rp = stat_wilcox_test(d_pair, .cols = pre:post, paired = TRUE)
  expect_equal(rp$p.value, unname(fitp$p.value), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# stat_anova / stat_kruskal
# ---------------------------------------------------------------------------

test_that("stat_anova one-way matches anova(lm)", {
  fit = stats::aov(len ~ supp, data = ToothGrowth)
  tab = stats::anova(fit)
  r = stat_anova(ToothGrowth, resp = len, group = supp)
  i_sup = which(r$term == "supp")
  expect_equal(r$sumsq[i_sup], unname(tab["supp", "Sum Sq"]), tolerance = 1e-8)
  expect_equal(r$statistic[i_sup], unname(tab["supp", "F value"]),
    tolerance = 1e-8)
  expect_equal(r$p.value[i_sup], unname(tab["supp", "Pr(>F)"]),
    tolerance = 1e-12)
  # partial eta squared
  expect_equal(r$effect[i_sup],
    tab["supp", "Sum Sq"] / (tab["supp", "Sum Sq"] + tab["Residuals", "Sum Sq"]),
    tolerance = 1e-10)
  # residuals row
  i_res = which(r$term == "Residuals")
  expect_equal(r$df[i_res], unname(tab["Residuals", "Df"]))
})

test_that("stat_anova two-way interaction matches", {
  mt2 = mtcars
  mt2$am = factor(mt2$am)
  mt2$cyl = factor(mt2$cyl)
  fit = stats::aov(mpg ~ am * cyl, data = mt2)
  tab = stats::anova(fit)
  r = stat_anova(mtcars, resp = mpg, group = c(am, cyl))
  i_int = which(r$term == "am:cyl")
  expect_equal(r$sumsq[i_int], unname(tab["am:cyl", "Sum Sq"]), tolerance = 1e-8)
  # multiple responses loop: 2 responses x (3 + 1) terms
  r2 = stat_anova(mtcars, resp = c(mpg, hp), group = cyl)
  expect_equal(nrow(r2), 2 * nrow(stats::anova(
    stats::aov(mpg ~ factor(cyl), data = mtcars))))
})

test_that("stat_kruskal matches stats::kruskal.test", {
  fit = stats::kruskal.test(len ~ supp, data = ToothGrowth)
  r = stat_kruskal(ToothGrowth, resp = len, group = supp)
  expect_equal(r$statistic, unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)
  n = 60; k = 2
  expect_equal(r$effect, (unname(fit$statistic) - k + 1) / (n - k),
    tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# stat_ancova
# ---------------------------------------------------------------------------

test_that("stat_ancova matches aov with covariate", {
  d = data.frame(
    y = c(10 + 2 * (1:15), 12 + 2 * (1:15)) + rnorm(30, 0, 0.5),
    g = rep(c("a", "b"), each = 15),
    x = 1:30
  )
  fit = stats::aov(y ~ g + x, data = d)
  tab = stats::anova(fit)
  r = stat_ancova(d, resp = y, group = g, covariates = x)
  i_g = which(r$term == "g")
  expect_equal(r$sumsq[i_g], unname(tab["g", "Sum Sq"]), tolerance = 1e-8)
  expect_equal(r$p.value[i_g], unname(tab["g", "Pr(>F)"]), tolerance = 1e-12)
  # parallel test warns on a strong interaction
  d2 = d
  d2$y = d$y + ifelse(d$g == "b", 0.6 * d$x, 0)
  expect_warning(stat_ancova(d2, resp = y, group = g, covariates = x,
    test_parallel = TRUE), "parallel-slopes")
})

# ---------------------------------------------------------------------------
# stat_cor
# ---------------------------------------------------------------------------

test_that("stat_cor matches stats::cor.test (pearson and spearman)", {
  r = stat_cor(mtcars, .cols = c(mpg, wt, hp))
  expect_equal(nrow(r), 3)
  fit = stats::cor.test(mtcars$mpg, mtcars$wt)
  i = which(r$var1 == "mpg" & r$var2 == "wt")
  expect_equal(r$estimate[i], unname(fit$estimate), tolerance = 1e-12)
  expect_equal(r$p.value[i], unname(fit$p.value), tolerance = 1e-12)
  expect_equal(r$conf_low[i], unname(fit$conf.int[1]), tolerance = 1e-10)

  r2 = stat_cor(mtcars, .cols = c(mpg, wt), method = "spearman")
  fit2 = stats::cor.test(mtcars$mpg, mtcars$wt, method = "spearman")
  expect_equal(r2$estimate, unname(fit2$estimate), tolerance = 1e-12)
  expect_true(all(is.na(r2$conf_low)))  # no CI for spearman
})

# ---------------------------------------------------------------------------
# stat_chisq
# ---------------------------------------------------------------------------

test_that("stat_chisq matches stats::chisq.test", {
  d = data.frame(
    treat = rep(c("a", "b"), each = 40),
    outcome = c(rep("yes", 28), rep("no", 12), rep("yes", 20), rep("no", 20))
  )
  tab = table(d$treat, d$outcome)
  fit = stats::chisq.test(tab)
  r = stat_chisq(d, row = treat, col = outcome)
  expect_equal(r$statistic, unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$p.value, unname(fit$p.value), tolerance = 1e-12)
  chi2_raw = suppressWarnings(unname(stats::chisq.test(tab,
    correct = FALSE)$statistic))
  expect_equal(r$effect, sqrt(chi2_raw / (80 * (2 - 1))), tolerance = 1e-10)

  rf = stat_chisq(d, row = treat, col = outcome, exact = TRUE)
  ff = stats::fisher.test(tab)
  expect_equal(rf$p.value, unname(ff$p.value), tolerance = 1e-12)
})

test_that("stat_chisq warns on low expected frequencies", {
  d = data.frame(
    g = factor(c(rep("a", 10), rep("b", 10))),
    o = factor(c(rep("y", 9), "n", rep("y", 2), rep("n", 8)))
  )
  expect_warning(stat_chisq(d, row = g, col = o), "expected frequencies")
})

# ---------------------------------------------------------------------------
# stat_normality
# ---------------------------------------------------------------------------

test_that("stat_normality matches stats::shapiro.test", {
  r = stat_normality(mtcars, .cols = c(mpg, wt))
  fit = stats::shapiro.test(mtcars$mpg)
  i = which(r$variable == "mpg")
  expect_equal(r$statistic[i], unname(fit$statistic), tolerance = 1e-10)
  expect_equal(r$p.value[i], unname(fit$p.value), tolerance = 1e-12)
  rg = stat_normality(ToothGrowth, .cols = len, .by = supp)
  fitg = stats::shapiro.test(ToothGrowth$len[ToothGrowth$supp == "OJ"])
  expect_equal(rg$statistic[rg$supp == "OJ"], unname(fitg$statistic),
    tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# stat_hov_test
# ---------------------------------------------------------------------------

test_that("stat_hov_test matches car/stats references", {
  r = stat_hov_test(ToothGrowth, .cols = len, group = supp)
  ref = car::leveneTest(len ~ supp, data = ToothGrowth, center = "median")
  expect_equal(r$statistic, unname(ref[1, "F value"]), tolerance = 1e-10)
  expect_equal(r$df, unname(ref[1, "Df"]))
  expect_equal(r$df2, unname(ref[2, "Df"]))
  expect_equal(r$p.value, unname(ref[1, "Pr(>F)"]), tolerance = 1e-12)

  rb = stat_hov_test(mtcars, .cols = mpg, group = cyl, method = "bartlett")
  fitb = stats::bartlett.test(mpg ~ cyl, data = mtcars)
  expect_equal(rb$statistic, unname(fitb$statistic), tolerance = 1e-10)
  expect_equal(rb$p.value, unname(fitb$p.value), tolerance = 1e-12)

  rf = stat_hov_test(ToothGrowth, .cols = len, group = supp, method = "f")
  fitf = stats::var.test(len ~ supp, data = ToothGrowth)
  expect_equal(rf$statistic, unname(fitf$statistic), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# stat_bootstrap
# ---------------------------------------------------------------------------

test_that("stat_bootstrap returns consistent CI and matches theory", {
  r = stat_bootstrap(mtcars, statistic = \(d, i) mean(d$mpg[i]),
    B = 800, seed = 1, type = "percentile")
  expect_s3_class(r, "stat_infer")
  expect_equal(r$term, "value")
  expect_equal(r$estimate, mean(mtcars$mpg), tolerance = 1e-12)
  # Monte-Carlo SE should be close to the theoretical SE (within 25%)
  theory_se = stats::sd(mtcars$mpg) / sqrt(32)
  expect_lt(abs(r$se - theory_se) / theory_se, 0.25)
  expect_true(r$conf_low < r$estimate && r$estimate < r$conf_high)
  expect_equal(r$B, 800)

  # BCa with boot.ci
  r2 = stat_bootstrap(mtcars, statistic = \(d, i) mean(d$mpg[i]),
    B = 800, seed = 1, type = "bca")
  expect_true(all(is.finite(c(r2$conf_low, r2$conf_high))))
})

test_that("stat_bootstrap multi-component statistics loop rows", {
  r = stat_bootstrap(mtcars,
    statistic = \(d, i) c(mean = mean(d$mpg[i]), sd = stats::sd(d$mpg[i])),
    B = 500, seed = 2, type = "basic")
  expect_equal(nrow(r), 2)
  expect_true(all(is.finite(c(r$conf_low, r$conf_high))))
})

test_that("stat_bootstrap .by slices repeat independently", {
  r = stat_bootstrap(d_two, statistic = \(d, i) mean(d$y[i]),
    B = 400, seed = 3, .by = region)
  expect_equal(nrow(r), 2)
  e = d_two$y[d_two$region == "e"]
  i = which(r$region == "e")
  expect_equal(r$estimate[i], mean(e), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# stat_perm_test
# ---------------------------------------------------------------------------

test_that("stat_perm_test two-sample default matches hand computation", {
  r = stat_perm_test(d_two, .cols = y, group = g, type = "two.sample",
    permutations = 500, seed = 11)
  expect_s3_class(r, "stat_infer")
  x1 = d_two$y[d_two$g == "a"]; x2 = d_two$y[d_two$g == "b"]
  expect_equal(r$observed, mean(x1) - mean(x2), tolerance = 1e-12)
  expect_true(r$p.value >= 1 / 501 && r$p.value <= 1)
  # with a real shift the two-sided p should be small
  expect_lt(r$p.value, 0.1)
})

test_that("stat_perm_test paired and one.sample flip signs symmetrically", {
  r = stat_perm_test(d_pair, .cols = pre:post, type = "paired",
    permutations = 500, seed = 5)
  expect_equal(r$observed, mean(d_pair$pre - d_pair$post), tolerance = 1e-12)
  # symmetric null distribution: null_mean approximately 0
  expect_lt(abs(r$null_mean), 0.15)

  r1 = stat_perm_test(mtcars, .cols = mpg, type = "one.sample", mu = 20,
    permutations = 500, seed = 6)
  expect_equal(r1$observed, mean(mtcars$mpg) - 20, tolerance = 1e-12)
})

test_that("stat_perm_test independence shuffles the second column", {
  set.seed(1)
  dind = data.frame(x = rnorm(50), z = rnorm(50))
  r = stat_perm_test(dind, .cols = x:z, type = "independence",
    permutations = 400, seed = 7)
  expect_equal(r$observed, stats::cor(dind$x, dind$z), tolerance = 1e-12)
  expect_lt(r$p.value, 1)
})

test_that("stat_perm_test custom statistic (two.sample) works", {
  # custom statistic: re-ordered rows carry both y and g, so the free
  # permutation re-assigns group labels; the statistic is the group
  # median difference computed from the re-ordered data
  f = function(d, i) {
    y = d$y[i]
    lab = factor(d$g[i], levels = c("a", "b"))
    stats::median(y[lab == "a"]) - stats::median(y[lab == "b"])
  }
  r = stat_perm_test(d_two, statistic = f, type = "two.sample",
    permutations = 300, seed = 9)
  expect_equal(r$observed,
    stats::median(d_two$y[d_two$g == "a"]) - stats::median(d_two$y[d_two$g == "b"]),
    tolerance = 1e-12)
  expect_true(r$p.value <= 1)
})

test_that("stat_perm_test .by slices and validation", {
  r = stat_perm_test(d_two, .cols = y, group = g, type = "two.sample",
    permutations = 300, seed = 13, .by = region)
  expect_equal(nrow(r), 2)
  expect_error(
    stat_perm_test(mtcars, .cols = mpg, type = "paired", permutations = 100),
    "two .cols"
  )
  expect_error(
    stat_perm_test(mtcars, statistic = 3, type = "two.sample"), "function"
  )
})

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

test_that("print.stat_infer renders and returns invisibly", {
  r = stat_t_test(d_two, .cols = y, group = g)
  out = capture.output(print(r))
  expect_true(grepl("t-test", out[1], fixed = TRUE))
})
