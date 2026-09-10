# Follow-up targeted checks for the audit report
pkgload::load_all("E:/Github/mathmodels", quiet = TRUE)

cat("=== A. stat_ancova test_parallel: ONLY first response checked? ===\n")
set.seed(5)
n = 16
dpar = data.frame(
  # response 1: strictly parallel slopes (no interaction)
  y1 = c(10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 18, 19) + rnorm(n, 0, .3),
  # response 2: strongly NON-parallel (group b slope +1.2)
  y2 = c(10, 11, 12, 13, 14, 15, 16, 17,
         10, 12.2, 14.4, 16.6, 18.8, 21, 23.2, 25.4) + rnorm(n, 0, .3),
  g = factor(rep(c("a", "b"), each = 8)),
  x = c(1:8, 1:8)
)
res = tryCatch({
  withCallingHandlers(
    stat_ancova(dpar, resp = c(y1, y2), group = g, covariates = x, test_parallel = TRUE),
    warning = function(w) { cat("  WARNING fired:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  "no-warning"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  resp = c(y1, y2) [only y2 non-parallel] ->", res, "\n")
res2 = tryCatch({
  withCallingHandlers(
    stat_ancova(dpar, resp = y2, group = g, covariates = x, test_parallel = TRUE),
    warning = function(w) { cat("  WARNING fired:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  "no-warning"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  resp = y2 alone ->", res2, "\n")
cat("  ==> if first line = no-warning and second warns: test_parallel skips responses beyond the first (L655)\n\n")

cat("=== B. stat_perm_test custom independence + .by: silent degenerate p ===\n")
dby = data.frame(x = 1:8, z = c(1, 3, 2, 5, 4, 7, 6, 8) + 0.5, s = rep(c("p", "q"), 4))
f_ind = function(d, i) cor(d$x, d$z)
r = stat_perm_test(dby, statistic = f_ind, type = "independence",
  permutations = 500, seed = 1, .by = s)
print(as.data.frame(r))
cat("  cor(x, z) overall =", signif(cor(dby$x, dby$z), 4),
    "; shuffled column is s (last), statistic invariant -> null_sd ~ 0, p ~ 1\n\n")

cat("=== C. stat_perm_test custom two.sample: package's own doc/test example trap ===\n")
dperm = data.frame(y = c(3, 5, 7, 9, 11, 2, 4, 6, 8, 10),
                   g = factor(rep(c("a", "b"), each = 5)))
# example copied from tests/testthat/test-stat_infer.R:373-378
f_test_example = function(d, i) {
  y = d$y[i]
  lab = factor(d$g[i], levels = c("a", "b"))
  stats::median(y[lab == "a"]) - stats::median(y[lab == "b"])
}
r2 = stat_perm_test(dperm, statistic = f_test_example, type = "two.sample",
  permutations = 2000, seed = 4)
cat(sprintf("  observed=%.4f  null_mean=%.6f  null_sd=%.6f  p=%.4f\n",
  r2$observed, r2$null_mean, r2$null_sd, r2$p.value))
cat("  ==> statistic invariant under permutation (indices applied to BOTH y and g):\n")
cat("      null degenerates to a point mass; p-value is meaningless (= 1)\n\n")

cat("=== D. stat_describe true n=1 and n=0 ===\n")
r1 = stat_describe(data.frame(x = 5), .cols = x)
cat("  n=1: n_used =", r1$n_used, " mean =", r1$mean, " sd =", r1$sd,
    " ci_low =", r1$ci_low, " skew =", r1$skew, "\n")
r0 = stat_describe(data.frame(x = numeric(0)), .cols = x)
cat("  n=0: n_used =", r0$n_used, " mean =", r0$mean, " (no crash)\n\n")

cat("=== E. t_test estimate direction (repro numbers for report) ===\n")
dt = data.frame(
  y = c(rnorm(15, 10, 2), rnorm(20, 11, 3)),
  g = factor(rep(c("a", "b"), c(15, 20))))
r = stat_t_test(dt, .cols = y, group = g)
cat(sprintf("  mean(a)=%.4f mean(b)=%.4f | statistic=%.4f (a-b) | estimate=%.4f (b-a) | effect=%.4f (a-b) | CI=[%.4f, %.4f] (a-b)\n",
  mean(dt$y[dt$g == "a"]), mean(dt$y[dt$g == "b"]),
  r$statistic, r$estimate, r$effect, r$conf_low, r$conf_high))

cat("\n=== F. stat_cor listwise vs pairwise-complete (repro) ===\n")
d = data.frame(a = 1:10, b = 1:10 + rnorm(10), cc = c(rep(1, 9), NA))
d$b[1] = NA
rc = stat_cor(d, .cols = a:cc)
print(as.data.frame(rc[, c("var1", "var2", "n_used")]))
cat("  pair a-b: complete rows = 9 (only row 1 has b=NA); reported n_used = 8 (row 10 cc=NA also dropped)\n")
cat("  stats::cor.test(a, b) n =", sum(!is.na(d$a) & !is.na(d$b)), "(pairwise-complete convention)\n")
