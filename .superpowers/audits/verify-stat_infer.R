# Verification probe: fix commit c8d4686 vs stat_infer audit findings F1-F10, F12
# Read-only on production code. Run:
#   D:\R-4.6.1\bin\Rscript.exe .superpowers/audits/verify-stat_infer.R

devtools::load_all("E:/Github/mathmodels", quiet = TRUE)

n_pass = 0L; n_fail = 0L; fails = character(0)
check = function(id, ok, detail = "") {
  ok = isTRUE(ok)
  if (ok) { n_pass <<- n_pass + 1L } else { n_fail <<- n_fail + 1L; fails <<- c(fails, id) }
  cat(if (ok) "PASS  " else "FAIL  ", id, if (nzchar(detail)) paste0("  [", detail, "]"), "\n", sep = "")
}
near = function(a, b, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
}
run_collect = function(expr) {
  ws = character(0)
  v = tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  list(value = if (inherits(v, "error")) NULL else v,
       warnings = ws,
       error = if (inherits(v, "error")) conditionMessage(v) else NULL)
}

cat("============================================================\n")
cat("F5: custom two.sample statistic (audit probe followup C)\n")
cat("============================================================\n")
dperm = data.frame(y = c(3, 5, 7, 9, 11, 2, 4, 6, 8, 10),
                   g = factor(rep(c("a", "b"), each = 5)))
# audit's EXACT statistic (old documented idiom, indexes by i)
f_audit = function(d, i) {
  y = d$y[i]; lab = factor(d$g[i], levels = c("a", "b"))
  stats::median(y[lab == "a"]) - stats::median(y[lab == "b"])
}
set.seed(4)
r5 = run_collect(stat_perm_test(dperm, statistic = f_audit, group = g,
  type = "two.sample", permutations = 2000, seed = 4))
cat(sprintf("  observed=%.4f null_mean=%.4f null_sd=%.4f p=%.4f (audit old: null_sd=0, p=1)\n",
  r5$value$observed, r5$value$null_mean, r5$value$null_sd, r5$value$p.value))
check("F5.1 meaningful null: null_sd > 0", r5$value$null_sd > 0,
  sprintf("null_sd=%.6f", r5$value$null_sd))
check("F5.2 p < 1 (not degenerate)", r5$value$p.value < 1,
  sprintf("p=%.4f", r5$value$p.value))
check("F5.3 no invariance warning on the good statistic", length(r5$warnings) == 0)
check("F5.4 observed = median contrast a-b (=1)", near(r5$value$observed, 1))
# audit's original call omitted group -> must now STOP cleanly (was: silent p=1)
r5b = run_collect(stat_perm_test(dperm, statistic = f_audit,
  type = "two.sample", permutations = 200, seed = 4))
check("F5.5 no-group call now stops with a clear message",
  !is.null(r5b$error) && grepl("group", r5b$error, ignore.case = TRUE),
  paste0("error: ", r5b$error))
# invariant statistic -> warning instead of silent p=1
set.seed(10)
r5c = run_collect(stat_perm_test(dperm, statistic = function(d, i) stats::median(d$y),
  group = g, type = "two.sample", permutations = 200, seed = 10))
check("F5.6 label-independent statistic warns 'invariant'",
  length(r5c$warnings) == 1 && grepl("invariant", r5c$warnings[1]),
  paste0("warnings: ", if (length(r5c$warnings)) r5c$warnings[1] else "(none)"))

cat("\n============================================================\n")
cat("F6: custom independence — second .cols column (audit case 2 + case 1)\n")
cat("============================================================\n")
dby = data.frame(x = 1:8, z = c(1, 3, 2, 5, 4, 7, 6, 8) + 0.5,
                 s = rep(c("p", "q"), 4))
f_ind = function(d, i) stats::cor(d$x, d$z)
cat(sprintf("  cor(x, z) overall = %.4f (audit: 0.93)\n", stats::cor(dby$x, dby$z)))
# audit case 2, now with .cols naming the test columns: shuffle z within slices
set.seed(1)
r6 = run_collect(stat_perm_test(dby, .cols = c(x, z), statistic = f_ind,
  type = "independence", permutations = 500, seed = 1, .by = s))
print(as.data.frame(r6$value))
check("F6.1 .by slices: null varies (z shuffled; slice col untouched)",
  all(r6$value$null_sd > 0) && all(r6$value$p.value < 1),
  sprintf("null_sd = %.4f, %.4f", r6$value$null_sd[1], r6$value$null_sd[2]))
check("F6.2 no invariance warning", length(r6$warnings) == 0)
cp = stats::cor(dby$x[dby$s == "p"], dby$z[dby$s == "p"])
cq = stats::cor(dby$x[dby$s == "q"], dby$z[dby$s == "q"])
check("F6.3 per-slice observed = within-slice cor(x, z)",
  near(r6$value$observed[r6$value$s == "p"], cp) &&
  near(r6$value$observed[r6$value$s == "q"], cq))
# audit case 1 flavor: trailing junk column must not become the shuffle target
set.seed(42)
dci = data.frame(x = 1:8, junk = stats::rnorm(8),
                 z = c(1, 3, 2, 5, 4, 7, 6, 8) + 0.5)
set.seed(6)
r6b = run_collect(stat_perm_test(dci, .cols = c(x, z, junk), statistic = f_ind,
  type = "independence", permutations = 2000, seed = 6))
check("F6.4 .cols = (x, z, junk): shuffles 2nd col z -> null_sd > 0, p < 1",
  r6b$value$null_sd > 0 && r6b$value$p.value < 1,
  sprintf("null_sd=%.4f p=%.4f", r6b$value$null_sd, r6b$value$p.value))
set.seed(6)
r6c = run_collect(stat_perm_test(dci, .cols = c(x, junk, z), statistic = f_ind,
  type = "independence", permutations = 500, seed = 6))
check("F6.5 statistic ignoring the 2nd .cols col -> invariant warning",
  length(r6c$warnings) == 1 && grepl("invariant", r6c$warnings[1]))
# clear stops
r6d = run_collect(stat_perm_test(dby, .cols = x, statistic = f_ind,
  type = "independence", permutations = 100, seed = 1))
check("F6.6 single .cols column -> clear stop",
  !is.null(r6d$error) && grepl("two", r6d$error), paste0("error: ", r6d$error))
r6e = run_collect(stat_perm_test(dby, statistic = f_ind,
  type = "independence", permutations = 100, seed = 1, .by = s))
check("F6.7 audit's old no-.cols call -> clear stop (was silent p = 1)",
  !is.null(r6e$error) && grepl(".cols", r6e$error), paste0("error: ", r6e$error))
r6f = run_collect(stat_perm_test(dby, .cols = c(x, s), statistic = f_ind,
  type = "independence", permutations = 100, seed = 1, .by = s))
check("F6.8 2nd .cols col is the .by slice col -> clear stop",
  !is.null(r6f$error) && grepl("slice", r6f$error), paste0("error: ", r6f$error))

cat("\n============================================================\n")
cat("F1: stat_cor pairwise-complete rows (audit repro)\n")
cat("============================================================\n")
set.seed(202)
dF = data.frame(a = 1:10, b = 1:10 + rnorm(10), cc = c(rep(1, 9), NA))
dF$b[1] = NA
rc = stat_cor(dF, .cols = a:cc)
print(as.data.frame(rc[, c("var1", "var2", "n_used", "estimate", "p.value")]))
i_ab = which(rc$var1 == "a" & rc$var2 == "b")
i_ac = which(rc$var1 == "a" & rc$var2 == "cc")
i_bc = which(rc$var1 == "b" & rc$var2 == "cc")
check("F1.1 (a,b) n_used = 9 (was 8)", rc$n_used[i_ab] == 9)
check("F1.2 (a,cc) n_used = 9", rc$n_used[i_ac] == 9)
check("F1.3 (b,cc) n_used = 8", rc$n_used[i_bc] == 8)
fit_ab = stats::cor.test(dF$a, dF$b)
check("F1.4 (a,b) matches cor.test (estimate)",
  near(rc$estimate[i_ab], unname(fit_ab$estimate), 1e-12))
check("F1.5 (a,b) matches cor.test (p)", near(rc$p.value[i_ab], unname(fit_ab$p.value), 1e-12))
fit_bc = stats::cor.test(dF$b, dF$cc)
check("F1.6 (b,cc) matches cor.test (estimate)",
  near(rc$estimate[i_bc], unname(fit_bc$estimate), 1e-12))

cat("\n============================================================\n")
cat("F7: stat_describe per-variable deletion + own-NA n_missing\n")
cat("============================================================\n")
set.seed(1)
dd = data.frame(x = c(rnorm(50), NA, 3), z = rnorm(52)); dd$z[1] = NA
r_alone = stat_describe(dd, .cols = x)
r_both = stat_describe(dd, .cols = x:z)
cat(sprintf("  mean(x) alone = %.5f | x-row in describe(x:z) = %.5f (audit old: differed)\n",
  r_alone$mean[1], r_both$mean[r_both$variable == "x"]))
check("F7.1 describe(x) mean == describe(x:z) x-mean",
  near(r_alone$mean[1], r_both$mean[r_both$variable == "x"]))
check("F7.2 x n_used = 51 in both calls",
  r_alone$n_used[1] == 51 && r_both$n_used[r_both$variable == "x"] == 51)
zrow = r_both[r_both$variable == "z", ]
check("F7.3 z n_used = 51, n_missing = 1 (own NA; audit old: 2)",
  zrow$n_used == 51 && zrow$n_missing == 1)
dg = data.frame(x = c(1:4, NA, 6:10), g = rep(c("a", "b"), each = 5))
rg = stat_describe(dg, .cols = x, .by = g)
print(as.data.frame(rg[, c("g", "n_used", "n_missing", "mean")]))
check("F7.4 slice b (0 own NAs): n_missing = 0 (audit old: 5), n_used = 5",
  rg$n_missing[rg$g == "b"] == 0 && rg$n_used[rg$g == "b"] == 5)
check("F7.5 slice a: n_missing = 1, n_used = 4",
  rg$n_missing[rg$g == "a"] == 1 && rg$n_used[rg$g == "a"] == 4)

cat("\n============================================================\n")
cat("F3: stat_ancova Type II SS on unbalanced data (audit section 6)\n")
cat("============================================================\n")
set.seed(5)
dc = data.frame(
  y  = c(10, 11, 12, 13, 14, 15, 16, 17, 9, 11, 13, 15, 17, 19, 21, 23) + rnorm(16, 0, .4),
  g  = factor(rep(c("a", "b"), each = 8)),
  x  = c(1:8, 1:8) + rnorm(16, 0, .2),
  x2 = rnorm(16)
)
dc2 = dc[1:14, ]
fit2 = stats::aov(y ~ g + x, data = dc2)
tI = stats::anova(fit2); tII = car::Anova(fit2, type = 2)
r3 = stat_ancova(dc2, resp = y, group = g, covariates = x)
ss_ours = r3$sumsq[r3$term == "g"]; p_ours = r3$p.value[r3$term == "g"]
cat(sprintf("  mathmodels SS_g = %.4f (p = %.4f) | TypeI = %.4f (p = %.4f) | car TypeII = %.4f (p = %.4f)\n",
  ss_ours, p_ours, tI["g", "Sum Sq"], tI["g", "Pr(>F)"],
  tII["g", "Sum Sq"], tII["g", "Pr(>F)"]))
cat(sprintf("  (audit reported: old mathmodels TypeI = 0.7511 vs TypeII = 8.4880)\n"))
check("F3.1 group SS matches car::Anova(aov, type = 2)",
  near(ss_ours, unname(tII["g", "Sum Sq"])))
check("F3.2 group p matches car type = 2", near(p_ours, unname(tII["g", "Pr(>F)"]), 1e-10))
check("F3.3 semantic change: differs from Type I on this unbalanced data",
  !near(ss_ours, unname(tI["g", "Sum Sq"]), 1e-6))
check("F3.4 internal: meansq = sumsq/df for all rows",
  all(near(r3$meansq, r3$sumsq / r3$df)))
ssr = r3$sumsq[r3$term == "Residuals"]
nr = which(r3$term != "Residuals")
check("F3.5 internal: partial eta^2 = ss/(ss + ss_res)",
  all(near(r3$effect[nr], r3$sumsq[nr] / (r3$sumsq[nr] + ssr))))
r3b = stat_ancova(dc, resp = y, group = g, covariates = x)
tIIb = car::Anova(stats::aov(y ~ g + x, data = dc), type = 2)
check("F3.6 balanced case matches car type 2",
  near(r3b$sumsq[r3b$term == "g"], unname(tIIb["g", "Sum Sq"])))
r3c = stat_ancova(dc, resp = y, group = g, covariates = c(x, x2))
tIIc = car::Anova(stats::aov(y ~ g + x + x2, data = dc), type = 2)
check("F3.7 two-covariate case matches car type 2",
  near(r3c$sumsq[r3c$term == "x"], unname(tIIc["x", "Sum Sq"])) &&
  near(r3c$sumsq[r3c$term == "x2"], unname(tIIc["x2", "Sum Sq"])))

cat("\n============================================================\n")
cat("F2: stat_t_test estimate direction = group1 - group2\n")
cat("============================================================\n")
set.seed(2)
dt = data.frame(y = c(rnorm(15, 10, 2), rnorm(20, 11, 3)),
                g = factor(rep(c("a", "b"), c(15, 20))))
r2_ = stat_t_test(dt, .cols = y, group = g)
ma = mean(dt$y[dt$g == "a"]); mb = mean(dt$y[dt$g == "b"])
cat(sprintf("  mean(a)=%.4f mean(b)=%.4f | estimate=%.4f | statistic=%.4f | CI=[%.4f, %.4f]\n",
  ma, mb, r2_$estimate, r2_$statistic, r2_$conf_low, r2_$conf_high))
check("F2.1 estimate = mean(a) - mean(b)", near(r2_$estimate, ma - mb, 1e-10))
check("F2.2 estimate sign consistent with statistic",
  sign(r2_$estimate) == sign(r2_$statistic))
check("F2.3 estimate1/estimate2 = group means",
  near(r2_$estimate1, ma, 1e-10) && near(r2_$estimate2, mb, 1e-10))
check("F2.4 CI (a-b direction) brackets estimate",
  r2_$conf_low <= r2_$estimate && r2_$estimate <= r2_$conf_high)

cat("\n============================================================\n")
cat("F4: test_parallel checks every response (audit probe followup A)\n")
cat("============================================================\n")
set.seed(5)
n = 16
dpar = data.frame(
  y1 = c(10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17, 18, 19) + rnorm(n, 0, .3),
  y2 = c(10, 11, 12, 13, 14, 15, 16, 17,
         10, 12.2, 14.4, 16.6, 18.8, 21, 23.2, 25.4) + rnorm(n, 0, .3),
  g = factor(rep(c("a", "b"), each = 8)),
  x = c(1:8, 1:8)
)
r4 = run_collect(stat_ancova(dpar, resp = c(y1, y2), group = g, covariates = x,
  test_parallel = TRUE))
for (w in r4$warnings) cat("  warning:", w, "\n")
check("F4.1 multi-response call warns; only y2 non-parallel (audit old: silent)",
  length(r4$warnings) >= 1 && any(grepl("parallel-slopes", r4$warnings)))
check("F4.2 warning names the violating response y2", any(grepl("y2", r4$warnings)))

cat("\n============================================================\n")
cat("F8 / F10 / F12\n")
cat("============================================================\n")
# F8: one-sample wilcox effect
fit8 = stats::wilcox.test(mtcars$mpg, mu = 20)
r8 = stat_wilcox_test(mtcars, .cols = mpg, mu = 20)
dv8 = mtcars$mpg - 20; n08 = sum(dv8 != 0); s_tot = n08 * (n08 + 1) / 2
rbc_hand = (2 * unname(fit8$statistic) - s_tot) / s_tot
cat(sprintf("  one-sample rbc = %.4f (hand = %.4f); audit old: NA\n", r8$effect, rbc_hand))
check("F8.1 one-sample effect non-NA", !is.na(r8$effect))
check("F8.2 effect = signed-rank rbc vs mu (hand check)", near(r8$effect, rbc_hand, 1e-10))
check("F8.3 p still matches stats::wilcox.test",
  near(r8$p.value, unname(fit8$p.value), 1e-12))
r8b = run_collect(stat_wilcox_test(data.frame(y = rep(5, 10)), .cols = y, mu = 5))
check("F8.4 all-zero deviations -> NA effect, no crash",
  is.null(r8b$error) && is.na(r8b$value$effect),
  paste0("error: ", r8b$error))
# F10: f-test df2
rf = stat_hov_test(ToothGrowth, .cols = len, group = supp, method = "f")
fitf = stats::var.test(len ~ supp, data = ToothGrowth)
cat(sprintf("  f-test df = %s df2 = %s (audit old: df2 = NA); var.test params = %s %s\n",
  rf$df, rf$df2, unname(fitf$parameter[1]), unname(fitf$parameter[2])))
check("F10.1 df2 present (was NA)", !is.na(rf$df2))
check("F10.2 df/df2 match var.test parameters",
  near(rf$df, unname(fitf$parameter[1])) && near(rf$df2, unname(fitf$parameter[2])))
# F12: kendall ties -> one clean package warning
dk = data.frame(x = c(1, 2, 2, 3, 4, 5, 5, 6), y = c(2, 1, 4, 3, 6, 5, 7, 8))
r12 = run_collect(stat_cor(dk, .cols = x:y, method = "kendall"))
cat("  kendall warnings:", if (length(r12$warnings))
  paste(r12$warnings, collapse = " | ") else "(none)", "\n")
check("F12.1 exactly one warning", length(r12$warnings) == 1)
check("F12.2 package-style ties message (ties present / asymptotic)",
  grepl("ties present", r12$warnings[1]) && grepl("asymptotic", r12$warnings[1]))
check("F12.3 raw stats message not leaked",
  !any(grepl("cannot compute exact", r12$warnings, ignore.case = TRUE)))
fit12 = suppressWarnings(stats::cor.test(dk$x, dk$y, method = "kendall"))
check("F12.4 estimate/p still match cor.test",
  near(r12$value$estimate, unname(fit12$estimate), 1e-12) &&
  near(r12$value$p.value, unname(fit12$p.value), 1e-12))

cat("\n============================================================\n")
cat(sprintf("SUMMARY: %d PASS / %d FAIL\n", n_pass, n_fail))
if (n_fail > 0) cat("Failed checks:", paste(fails, collapse = ", "), "\n")
cat("============================================================\n")
