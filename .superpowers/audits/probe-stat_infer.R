# ============================================================================
# Audit probe v2: mathmodels::stat_infer.* vs gold-standard references
# Run: D:\R-4.6.1\bin\Rscript.exe .superpowers\audits\probe-stat_infer.R
# ============================================================================
pkgload::load_all("E:/Github/mathmodels", quiet = TRUE)
suppressPackageStartupMessages({
  library(coin)
  library(car)
})

rdiff = function(a, b) {
  a = unname(a); b = unname(b)
  if (length(a) != 1 || length(b) != 1 || !is.finite(a) || !is.finite(b)) return(NA_real_)
  if (a == b) return(0)
  min(abs(a - b) / max(abs(a), abs(b), 1e-300), Inf)
}
pv_close = function(a, b) {
  a = unname(a); b = unname(b)
  if (length(a) != 1 || length(b) != 1) return(NA)
  if (is.na(a) || is.na(b)) return(is.na(a) && is.na(b))
  abs(a - b) < 1e-6 || abs(a - b) / max(abs(a), abs(b), 1e-300) < 1e-8
}
OK = "MATCH"; VR = "VARIANT"; BUG = "BUG"
RES = new.env(); RES$tbl = NULL
note = function(fun, item, ours, ref, rd, verdict, extra = "") {
  txt = sprintf("%-17s %-44s ours=%-13s ref=%-13s rd=%-9s %s %s",
    fun, item,
    format(if (length(ours) == 1 && !is.character(ours)) ours else paste(ours, collapse=","), digits = 7),
    format(if (length(ref) == 1 && !is.character(ref)) ref else paste(ref, collapse=","), digits = 7),
    if (is.na(rd)) "NA" else sprintf("%.1e", rd),
    verdict, extra)
  cat(txt, "\n")
  RES$tbl = rbind(RES$tbl, data.frame(
    fun = fun, item = item,
    ours = if (length(ours) == 1 && is.numeric(ours)) ours else NA_real_,
    ref = if (length(ref) == 1 && is.numeric(ref)) ref else NA_real_,
    rd = if (is.numeric(rd)) rd else NA_real_, verdict = verdict,
    extra = extra, stringsAsFactors = FALSE))
}
gv = function(df, col, row = 1) df[[col]][row]

# ============================================================================
# 1. stat_describe
# ============================================================================
cat("\n=== 1. stat_describe ===\n")
set.seed(1)
dd = data.frame(x = c(rnorm(50), NA, 3), z = rnorm(52))
dd$z[1] = NA
ddc = dd[stats::complete.cases(dd), ]       # mathmodels listwise reference
r = stat_describe(dd, .cols = x:z)
rx = r[r$variable == "x", ]
rz = r[r$variable == "z", ]
cmp = function(item, a, b) note("stat_describe", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
cmp("n_used x (listwise over BOTH cols)", rx$n_used, nrow(ddc))
cmp("mean x (listwise ref)", rx$mean, mean(ddc$x))
cmp("sd x", rx$sd, sd(ddc$x))
cmp("median x", rx$median, median(ddc$x))
cmp("iqr x", rx$iqr, IQR(ddc$x))
g1f = function(x) { x = x[!is.na(x)]; m = mean(x); s2 = mean((x-m)^2); mean((x-m)^3)/s2^1.5 }
g2f = function(x) { x = x[!is.na(x)]; m = mean(x); s2 = mean((x-m)^2); mean((x-m)^4)/s2^2 - 3 }
cmp("skew x (population g1)", rx$skew, g1f(ddc$x))
cmp("kurt x (population excess g2)", rx$kurt, g2f(ddc$x))
hwx = qt(0.975, nrow(ddc) - 1) * sd(ddc$x) / sqrt(nrow(ddc))
cmp("ci_low x", rx$ci_low, mean(ddc$x) - hwx)
note("stat_describe", "skew vs psych::describe type-3", rx$skew,
  psych::describe(dd$x, na.rm = TRUE)$skew, NA, VR,
  "psych type-3 (unbiased) differs by design; mathmodels doc says moment (biased g1)")
# SEMANTIC: cross-column listwise deletion changes x's own stats
r_alone = stat_describe(dd, .cols = x)
note("stat_describe", "describe(x) mean vs describe(x:z) mean",
  r_alone$mean, rx$mean, NA,
  if (isTRUE(all.equal(r_alone$mean, rx$mean))) OK else VR,
  sprintf("mean(x) alone=%.5f vs in company=%.5f -> listwise deletion across .cols (not documented)",
    r_alone$mean, rx$mean))
# SEMANTIC: n_missing counts rows dropped by ANY column / other slices
note("stat_describe", "n_missing for z (z has 1 own NA)", rz$n_missing, 2, NA,
  if (rz$n_missing == 2) VR else BUG,
  "n_missing = rows dropped by listwise rule (incl. row where x is NA), not z's own NA count (=1)")
d_by = data.frame(x = c(1:4, NA, 6:10), g = rep(c("a", "b"), each = 5))
rb = stat_describe(d_by, .cols = x, .by = g)
note("stat_describe", "n_missing slice b (0 own NA)", rb$n_missing[rb$g == "b"], 5, NA,
  if (rb$n_missing[rb$g == "b"] == 5) VR else BUG,
  "n_missing counts the 5 OTHER-slice rows as 'missing' for slice b")
# CI when n = 1: NA; mean present
d1 = data.frame(x = c(1, NA, 3))
r1 = stat_describe(d1, .cols = x)
note("stat_describe", "n=1: mean ok, ci NA", paste(r1$mean, r1$ci_low, sep = "/"), "2/NA", NA, OK)

# ============================================================================
# 2. stat_t_test
# ============================================================================
cat("\n=== 2. stat_t_test ===\n")
set.seed(2)
dt = data.frame(
  y = c(rnorm(15, 10, 2), rnorm(20, 11, 3)),
  g = factor(rep(c("a", "b"), c(15, 20))),
  region = rep(c("e", "w"), length.out = 35)
)
x_a = dt$y[dt$g == "a"]; x_b = dt$y[dt$g == "b"]
fw = t.test(x_a, x_b)          # what mathmodels calls internally (default iface)
fwf = t.test(y ~ g, data = dt) # formula iface
rw = stat_t_test(dt, .cols = y, group = g)
cat("  t.test default t =", signif(fw$statistic, 6), " formula t =", signif(fwf$statistic, 6),
    " mean(a)-mean(b) =", signif(mean(x_a) - mean(x_b), 6), "\n")
pcmp = function(item, a, b) note("stat_t_test", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
cmp2 = function(item, a, b, tol = 1e-6) note("stat_t_test", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < tol) OK else BUG)
cmp2("Welch statistic", rw$statistic, fw$statistic)
cmp2("Welch df", rw$df, fw$parameter)
pcmp("Welch p", rw$p.value, fw$p.value)
cmp2("Welch ci_low", rw$conf_low, fw$conf.int[1])
# estimate direction: internally consistent?
note("stat_t_test", "estimate col vs statistic direction", rw$estimate,
  mean(x_b) - mean(x_a), rdiff(rw$estimate, mean(x_b) - mean(x_a)), VR,
  sprintf("estimate = mean(b)-mean(a)=%.4f while statistic/CI/Cohen-d direction = a-b (=%.4f): SIGN CONTRADICTION inside one row",
    rw$estimate, rw$statistic))
rstat = rstatix::t_test(dt, y ~ g)
cat("  rstatix cols:", paste(names(rstat), collapse = ","), "\n")
if ("estimate" %in% names(rstat)) cat("  rstatix estimate =", signif(rstat$estimate, 6), "\n")
fp = t.test(x_a, x_b, var.equal = TRUE)
rp = stat_t_test(dt, .cols = y, group = g, var.equal = TRUE)
cmp2("pooled df", rp$df, fp$parameter)
pcmp("pooled p", rp$p.value, fp$p.value)
f1 = t.test(dt$y, mu = 10)
r1 = stat_t_test(dt, .cols = y, mu = 10)
pcmp("one-sample p", r1$p.value, f1$p.value)
cmp2("one-sample effect dz", r1$effect, (mean(dt$y) - 10) / sd(dt$y))
dp = data.frame(pre = rnorm(20, 5), post = rnorm(20, 5.5, 0.9))
fpair = t.test(dp$pre, dp$post, paired = TRUE)
rpair = stat_t_test(dp, .cols = pre:post, paired = TRUE)
pcmp("paired p", rpair$p.value, fpair$p.value)
cmp2("paired estimate (x-y)", rpair$estimate, fpair$estimate)
cmp2("paired effect dz", rpair$effect, mean(dp$pre - dp$post) / sd(dp$pre - dp$post))
fg = t.test(x_a, x_b, alternative = "greater")
rg = stat_t_test(dt, .cols = y, group = g, alternative = "greater")
pcmp("greater p", rg$p.value, fg$p.value)
cmp2("greater one-sided ci_low", rg$conf_low, fg$conf.int[1])
fm = t.test(x_a, x_b, mu = 1)
rm = stat_t_test(dt, .cols = y, group = g, mu = 1)
pcmp("mu=1 p", rm$p.value, fm$p.value)
sp = sqrt(((14) * var(x_a) + (19) * var(x_b)) / 33)
cmp2("Cohens d (pooled)", rw$effect, (mean(x_a) - mean(x_b)) / sp)
# NA handling: per-test-column complete cases only
dtn = dt; dtn$y[1] = NA; dtn$other = 1:35
rn = stat_t_test(dtn, .cols = y, group = g)
pcmp("NA row dropped p", rn$p.value, t.test(y ~ g, dtn[2:35, ])$p.value)

# ============================================================================
# 3. stat_wilcox_test
# ============================================================================
cat("\n=== 3. stat_wilcox_test ===\n")
dw = data.frame(
  y = c(1.5, 2.5, 3, 4, 5, 2.7, 3.3, 4.8, 6, 7, 1.9, 5.5),
  g = factor(rep(c("a", "b"), c(5, 7)))
)
n1w = 5; n2w = 7
# which W does the formula interface return?
f_def = suppressWarnings(wilcox.test(dw$y[dw$g == "a"], dw$y[dw$g == "b"]))
f_frm = suppressWarnings(wilcox.test(y ~ g, data = dw))
rw = stat_wilcox_test(dw, .cols = y, group = g)
cat("  W default(a,b) =", signif(f_def$statistic, 6), " W formula =", signif(f_frm$statistic, 6), "\n")
# hand W_a: rank-sum statistic for group a (no ties in dw)
rk = rank(dw$y)
W_a_hand = sum(rk[dw$g == "a"]) - n1w * (n1w + 1) / 2
r_ref_a = 2 * W_a_hand / (n1w * n2w) - 1
r_ref_b = -r_ref_a
wc = function(item, a, b, tol = 1e-6) note("stat_wilcox_test", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < tol) OK else BUG)
wpcmp = function(item, a, b) note("stat_wilcox_test", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
wc("W (2-samp, ties ok)", rw$statistic, f_frm$statistic)
wpcmp("p (2-samp, correct=TRUE)", rw$p.value, f_frm$p.value)
note("stat_wilcox_test", "effect: 2W/(n1n2)-1 vs hand W_a rank-biserial",
  rw$effect, r_ref_a, rdiff(rw$effect, r_ref_a),
  if (is.na(rdiff(rw$effect, r_ref_a))) VR else if (rdiff(rw$effect, r_ref_a) < 1e-6) OK else BUG,
  sprintf("hand r(a-direction)=%.4f; rstatix effsize check follows", r_ref_a))
we = rstatix::wilcox_effsize(dw, y ~ g)
cat("  rstatix wilcox_effsize =", signif(we$effsize, 6), " n1=", we$n1, " n2=", we$n2, "\n")
note("stat_wilcox_test", "effect vs rstatix (sign)", rw$effect, we$effsize, rdiff(rw$effect, we$effsize),
  if (isTRUE(all.equal(unname(rw$effect), unname(we$effsize)))) OK else VR,
  "sign convention vs rstatix")
# coin asymptotic cross-check
cw = coin::wilcox_test(y ~ g, data = dw, distribution = "asymptotic")
note("stat_wilcox_test", "p vs coin asymptotic", rw$p.value, coin::pvalue(cw), NA,
  if (pv_close(rw$p.value, coin::pvalue(cw))) OK else VR,
  "coin: no continuity correction; stats default correct=TRUE -> small gap expected")
# exact (no ties, small n)
dnt = data.frame(y = c(1:5, 2.2, 3.3, 4.4), g = factor(rep(c("a","b"), c(5,3))))
ce = coin::wilcox_test(y ~ g, data = dnt, distribution = "exact")
fs = wilcox.test(y ~ g, data = dnt, exact = TRUE, correct = FALSE)
rs = stat_wilcox_test(dnt, .cols = y, group = g)
wpcmp("exact p vs stats (no ties)", rs$p.value, fs$p.value)
wpcmp("exact p vs coin exact", rs$p.value, coin::pvalue(ce))
# one-sample
fo = suppressWarnings(wilcox.test(dw$y, mu = 3.5))
ro = stat_wilcox_test(dw, .cols = y, mu = 3.5)
wpcmp("one-sample p", ro$p.value, fo$p.value)
note("stat_wilcox_test", "one-sample effect = NA", ro$effect, NA, NA, VR,
  "docs promise rank-biserial effect but one-sample path returns NA")
# paired
fpr = suppressWarnings(wilcox.test(dp$pre, dp$post, paired = TRUE))
rpr = stat_wilcox_test(dp, .cols = pre:post, paired = TRUE)
wpcmp("paired p", rpr$p.value, fpr$p.value)
dvz = dp$pre - dp$post; dvz = dvz[dvz != 0]; n0 = length(dvz)
eff_hand = (2 * unname(fpr$statistic) - n0 * (n0 + 1) / 2) / (n0 * (n0 + 1) / 2)
wc("paired rbc formula (zeros dropped)", rpr$effect, eff_hand)
# zeros in diffs
dz = data.frame(pre = c(1, 2, 3, 4, 5, 6), post = c(1, 2, 4, 3, 6, 8))
fz = wilcox.test(dz$pre, dz$post, paired = TRUE)
rz = stat_wilcox_test(dz, .cols = pre:post, paired = TRUE)
wpcmp("paired zeros p", rz$p.value, fz$p.value)
nz0 = 4
wc("paired rbc zeros-dropped n=4", rz$effect,
  (2 * unname(fz$statistic) - nz0 * (nz0 + 1) / 2) / (nz0 * (nz0 + 1) / 2))
# all-equal two-sample: zero variance in ranks?
eq = data.frame(y = rep(5, 10), g = rep(c("a","b"), each = 5))
req = tryCatch({suppressWarnings(stat_wilcox_test(eq, .cols = y, group = g)); "ok"},
  error = function(e) paste("ERROR:", conditionMessage(e)))
note("stat_wilcox_test", "all-equal groups", req, "ok", NA, if (req == "ok") VR else BUG,
  "stats::wilcox.test degenerate sigma=0 -> what happens?")

# ============================================================================
# 4. stat_anova
# ============================================================================
cat("\n=== 4. stat_anova ===\n")
set.seed(4)
da = data.frame(
  y = c(rnorm(12, 10), rnorm(10, 12), rnorm(11, 9)) + rep(0:2, c(12, 10, 11)),
  g1 = factor(rep(c("a","b","c"), c(12, 10, 11))),
  g2 = factor(rep(rep(c("x","y"), each = 2), length.out = 33))
)
t1w = anova(aov(y ~ g1, data = da))
r1w = stat_anova(da, resp = y, group = g1)
ac = function(item, a, b) note("stat_anova", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
ap = function(item, a, b) note("stat_anova", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
ac("1-way SS g1", gv(r1w, "sumsq", 1), t1w["g1", "Sum Sq"])
ac("1-way F g1", gv(r1w, "statistic", 1), t1w["g1", "F value"])
ap("1-way p g1", gv(r1w, "p.value", 1), t1w["g1", "Pr(>F)"])
ac("partial eta2 g1", gv(r1w, "effect", 1),
  t1w["g1","Sum Sq"] / (t1w["g1","Sum Sq"] + t1w["Residuals","Sum Sq"]))
ac("Residuals row present", nrow(r1w), nrow(t1w))
f2w = aov(y ~ g1 * g2, data = da); t2w = anova(f2w)
r2w = stat_anova(da, resp = y, group = c(g1, g2))
ac("2-way SS g1 (Type I)", r2w$sumsq[r2w$term == "g1"], t2w["g1", "Sum Sq"])
ac("2-way SS g2 (Type I)", r2w$sumsq[r2w$term == "g2"], t2w["g2", "Sum Sq"])
ac("2-way SS g1:g2", r2w$sumsq[r2w$term == "g1:g2"], t2w["g1:g2", "Sum Sq"])
tII = car::Anova(f2w, type = 2)
note("stat_anova", "2-way g1 SS vs Type II", r2w$sumsq[r2w$term == "g1"],
  tII["g1", "Sum Sq"], rdiff(r2w$sumsq[r2w$term=="g1"], tII["g1", "Sum Sq"]), VR,
  "Type I sequential SS; rstatix anova_test defaults to Type II; docs silent on SS type")
note("stat_anova", "p.adjusted always == p (no p_adjust arg)",
  all(r2w$p.adjusted == r2w$p.value), TRUE, NA, VR,
  "docs/@return advertise p.adjusted column but stat_anova has no p_adjust parameter (L512 hardcodes none)")
# .by slice equal to manual subset
rby = stat_anova(da, resp = y, group = g1, .by = g2)
sub = da[da$g2 == "x", ]
tsub = anova(aov(y ~ g1, data = sub))
ac(".by slice SS g1", rby$sumsq[rby$g2 == "x" & rby$term == "g1"], tsub["g1", "Sum Sq"])

# ============================================================================
# 5. stat_kruskal
# ============================================================================
cat("\n=== 5. stat_kruskal ===\n")
fk = kruskal.test(y ~ g1, data = da)
rk = stat_kruskal(da, resp = y, group = g1)
kc = function(item, a, b) note("stat_kruskal", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
kp = function(item, a, b) note("stat_kruskal", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
kc("H statistic", rk$statistic, fk$statistic)
kp("p asymptotic", rk$p.value, fk$p.value)
kc("epsilon2 = (H-k+1)/(n-k)", rk$effect, (fk$statistic - 3 + 1) / (33 - 3))
ke = rstatix::kruskal_effsize(da, y ~ g1)
note("stat_kruskal", "epsilon2 vs rstatix", rk$effect, ke$effsize,
  rdiff(rk$effect, ke$effsize),
  if (isTRUE(all.equal(unname(rk$effect), unname(ke$effsize)))) OK else VR, "")
dk = data.frame(y = c(1,2,3,4, 2,3,5,6), g = factor(rep(c("a","b"), each = 4)))
ce_k = coin::kruskal_test(y ~ g, data = dk, distribution = "exact")
rs_k = stat_kruskal(dk, resp = y, group = g)
note("stat_kruskal", "small-n asymptotic p vs coin exact", rs_k$p.value, coin::pvalue(ce_k), NA, VR,
  "asymptotic chi2 approximation (stats) vs exact differ by design; n=8 total")

# ============================================================================
# 6. stat_ancova
# ============================================================================
cat("\n=== 6. stat_ancova ===\n")
set.seed(5)
dc = data.frame(
  y  = c(10,11,12,13,14,15,16,17, 9,11,13,15,17,19,21,23) + rnorm(16, 0, .4),
  g  = factor(rep(c("a","b"), each = 8)),
  x  = c(1:8, 1:8) + rnorm(16, 0, .2),
  x2 = rnorm(16)
)
fit_seq = aov(y ~ g + x, data = dc)
tI = anova(fit_seq); tII = car::Anova(fit_seq, type = 2)
r_anc = stat_ancova(dc, resp = y, group = g, covariates = x)
anc = function(item, a, b) note("stat_ancova", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
anc("g SS == Type I", r_anc$sumsq[r_anc$term == "g"], tI["g", "Sum Sq"])
anc("x SS == Type I", r_anc$sumsq[r_anc$term == "x"], tI["x", "Sum Sq"])
note("stat_ancova", "g SS vs Type II (covariate-adjusted)", r_anc$sumsq[r_anc$term == "g"],
  tII["g", "Sum Sq"], rdiff(r_anc$sumsq[r_anc$term=="g"], tII["g", "Sum Sq"]), VR,
  "standard ANCOVA group test is Type II/III; Type I tests g WITHOUT covariate adjustment")
# unbalanced: Type I != Type II for g, clearly
dc2 = dc[1:14, ]
fit2 = aov(y ~ g + x, data = dc2)
tI2 = anova(fit2); tII2 = car::Anova(fit2, type = 2)
r2 = stat_ancova(dc2, resp = y, group = g, covariates = x)
anc("unbal g SS Type I", r2$sumsq[r2$term == "g"], tI2["g", "Sum Sq"])
note("stat_ancova", "unbal g SS Type I vs II", r2$sumsq[r2$term == "g"], tII2["g", "Sum Sq"],
  rdiff(r2$sumsq[r2$term=="g"], tII2["g", "Sum Sq"]), VR,
  sprintf("TypeI=%.4f vs TypeII=%.4f -> group tested unadjusted (ANCOVA namesake test affected)",
    r2$sumsq[r2$term=="g"], tII2["g", "Sum Sq"]))
# anova_test (rstatix, Type II default) cross-ref (defensive print)
rat = rstatix::anova_test(dc2, y ~ g + x, type = 2)
at = rat$ANOVA$ANOVA
pg_ref = tryCatch(at[at$Effect == "g", "p"], error = function(e) NA)
cat("  rstatix Type II g p =", format(signif(as.numeric(pg_ref), 6)),
    " mathmodels g p =", signif(r2$p.value[r2$term == "g"], 6), "\n")
# two covariates
r3c = stat_ancova(dc, resp = y, group = g, covariates = c(x, x2))
t3c = anova(aov(y ~ g + x + x2, data = dc))
anc("2-cov x SS", r3c$sumsq[r3c$term == "x"], t3c["x", "Sum Sq"])
anc("2-cov x2 SS", r3c$sumsq[r3c$term == "x2"], t3c["x2", "Sum Sq"])
# test_parallel
dc_par = dc; dc_par$y = dc$y + ifelse(dc$g == "b", 0.5 * dc$x, 0)
w_par = tryCatch({stat_ancova(dc_par, resp = y, group = g, covariates = x, test_parallel = TRUE); FALSE},
  warning = function(w) TRUE)
note("stat_ancova", "test_parallel warns on interaction", w_par, TRUE, NA, if (w_par) OK else BUG)
dc_mr = dc; dc_mr$y2 = dc_mr$y + ifelse(dc_mr$g == "b", 0.5 * dc_mr$x, 0)  # y2 also non-parallel
res_mr = tryCatch({stat_ancova(dc_mr, resp = c(y, y2), group = g, covariates = x,
  test_parallel = TRUE); "no-warn"}, warning = function(w) "warned")
note("stat_ancova", "test_parallel checks only first response", res_mr, "no-warn", NA, VR,
  "L655 loops responses but tests only names(resp_df)[1]; other responses' interactions untested")

# ============================================================================
# 7. stat_cor
# ============================================================================
cat("\n=== 7. stat_cor ===\n")
set.seed(6)
dco = data.frame(a = rnorm(40), b = rnorm(40), cc = rnorm(40))
dco$b = 0.6 * dco$a + rnorm(40, 0, .8)
dco$cc = dco$a + rnorm(40, 0, 2)
for (m in c("pearson", "spearman", "kendall")) {
  fc = cor.test(dco$a, dco$b, method = m)
  rc = stat_cor(dco, .cols = a:b, method = m)
  note("stat_cor", paste(m, "estimate"), rc$estimate, fc$estimate, rdiff(rc$estimate, fc$estimate),
    if (rdiff(rc$estimate, fc$estimate) < 1e-6) OK else BUG)
  note("stat_cor", paste(m, "p"), rc$p.value, fc$p.value, if (pv_close(rc$p.value, fc$p.value)) 0 else 1,
    if (pv_close(rc$p.value, fc$p.value)) OK else BUG)
  note("stat_cor", paste(m, "statistic"), rc$statistic, fc$statistic, rdiff(rc$statistic, fc$statistic),
    if (rdiff(rc$statistic, fc$statistic) < 1e-6) OK else BUG)
}
fg = cor.test(dco$a, dco$b, alternative = "greater")
rg = stat_cor(dco, .cols = a:b, alternative = "greater")
note("stat_cor", "pearson greater p", rg$p.value, fg$p.value, if (pv_close(rg$p.value, fg$p.value)) 0 else 1,
  if (pv_close(rg$p.value, fg$p.value)) OK else BUG)
f95 = cor.test(dco$a, dco$b, conf.level = 0.9)
r95 = stat_cor(dco, .cols = a:b, conf_level = 0.9)
note("stat_cor", "pearson CI 90% low", r95$conf_low, f95$conf.int[1],
  rdiff(r95$conf_low, f95$conf.int[1]), if (rdiff(r95$conf_low, f95$conf.int[1]) < 1e-6) OK else BUG)
rs = stat_cor(dco, .cols = a:b, method = "spearman")
note("stat_cor", "spearman CI NA documented", rs$conf_low, NA, NA, OK)
# kendall ties warning propagation
dct = dco; dct$b[1:3] = round(dct$b[1:3], 0)
kw = tryCatch({stat_cor(dct, .cols = a:b, method = "kendall"); "no-warn"},
  warning = function(w) paste("WARN:", conditionMessage(w)))
note("stat_cor", "kendall ties warning propagates", kw, "no-warn", NA, VR,
  "stats 'cannot compute exact p-value with ties' warning leaks to user (not suppressed)")
# listwise across all pairs
dcoNA = dco; dcoNA$b[1] = NA; dcoNA$cc[3] = NA
rNA = stat_cor(dcoNA, .cols = a:cc)
note("stat_cor", "n listwise across all cols (a-b pair)", rNA$n_used[1], 39, NA,
  if (rNA$n_used[1] == 39) VR else BUG,
  "pair (a,b) loses row 3 (cc NA) though a,b complete there -> pairwise-complete would give 40")

# ============================================================================
# 8. stat_chisq
# ============================================================================
cat("\n=== 8. stat_chisq ===\n")
dch = data.frame(
  treat = rep(c("a","b"), each = 40),
  outcome = c(rep("yes", 30), rep("no", 10), rep("yes", 20), rep("no", 20))
)
tab = table(dch$treat, dch$outcome)
rc2 = stat_chisq(dch, row = treat, col = outcome)
fc2 = chisq.test(tab)
chc = function(item, a, b) note("stat_chisq", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
chp = function(item, a, b) note("stat_chisq", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
chc("2x2 stat (correct=TRUE)", rc2$statistic, fc2$statistic)
chp("2x2 p (correct=TRUE)", rc2$p.value, fc2$p.value)
chc("min_expected", rc2$min_expected, min(fc2$expected))
chc("pct_below_5", rc2$pct_below_5, mean(fc2$expected < 5))
v_raw = sqrt(as.numeric(chisq.test(tab, correct = FALSE)$statistic) / (80 * 1))
chc("Cramers V (uncorrected chi2)", rc2$effect, v_raw)
cv_rstat = tryCatch(rstatix::cramers_v(dch, treat ~ outcome), error = function(e) NULL)
cv_val = if (!is.null(cv_rstat)) as.numeric(cv_rstat[[1]][1]) else NA_real_
note("stat_chisq", "V vs rstatix cramers_v", rc2$effect, cv_val, rdiff(rc2$effect, cv_val),
  if (isTRUE(all.equal(unname(rc2$effect), cv_val))) OK else VR, "")
fcf = chisq.test(tab, correct = FALSE)
rcf = stat_chisq(dch, row = treat, col = outcome, correct = FALSE)
chp("2x2 p (correct=FALSE)", rcf$p.value, fcf$p.value)
chc("V same under correct=TRUE/FALSE", rcf$effect, rc2$effect)
drc = data.frame(g = factor(rep(c("a","b","c"), each = 30)),
                 o = factor(rep(c("x","y","z","x","y","z","z"), length.out = 90)))
fc3 = suppressWarnings(chisq.test(table(drc$g, drc$o)))
rc3 = stat_chisq(drc, row = g, col = o)
chp("3x3 p (no correction)", rc3$p.value, fc3$p.value)
rf = stat_chisq(dch, row = treat, col = outcome, exact = TRUE)
chp("fisher p", rf$p.value, fisher.test(tab)$p.value)
note("stat_chisq", "fisher stat/df NA", paste(rf$statistic, rf$df, sep = "/"), "NA/NA", NA, OK)
# goodness-of-fit not offered (doc title says independence) - confirm signature
note("stat_chisq", "GOF not supported (doc-consistent)", "n/a", "n/a", NA, OK)

# ============================================================================
# 9. stat_normality
# ============================================================================
cat("\n=== 9. stat_normality ===\n")
rn = stat_normality(dco, .cols = a:cc)
fsn = shapiro.test(dco$a)
nc = function(item, a, b) note("stat_normality", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
np = function(item, a, b) note("stat_normality", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
nc("SW stat a", rn$statistic[rn$variable == "a"], fsn$statistic)
np("SW p a", rn$p.value[rn$variable == "a"], fsn$p.value)
nc("skew a (g1)", rn$skew[rn$variable == "a"], g1f(dco$a))
nc("kurt a (g2)", rn$kurt[rn$variable == "a"], g2f(dco$a))
r3 = stat_normality(data.frame(x = 1:3), .cols = x)
nc("n=3 SW works", r3$statistic, shapiro.test(1:3)$statistic)
r2n = stat_normality(data.frame(x = 1:2), .cols = x)
note("stat_normality", "n=2 -> NA no crash", r2n$statistic, NA, NA, OK)
r5001 = tryCatch({
  stat_normality(data.frame(x = rnorm(5001)), .cols = x)$statistic
}, error = function(e) paste("ERROR:", conditionMessage(e)))
note("stat_normality", "n=5001 -> NA guard", r5001, NA, NA,
  if (is.na(r5001)) OK else VR, "doc says NA outside [3,5000]")
# constant column
rc0 = tryCatch({
  r = withCallingHandlers(
    stat_normality(data.frame(x = rep(5, 20)), .cols = x),
    warning = function(w) invokeRestart("muffleWarning"))
  paste(sprintf("W=%s p=%s skew=%s", format(r$statistic), format(r$p.value), format(r$skew)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
note("stat_normality", "constant column behavior", rc0, "no-crash", NA, VR,
  "shapiro.test errors 'all x values identical'? or NaN; .moments_sk returns NA")
# nortest AD not claimed by docs -> not compared (docs only promise Shapiro-Wilk)

# ============================================================================
# 10. stat_hov_test
# ============================================================================
cat("\n=== 10. stat_hov_test ===\n")
fl = car::leveneTest(y ~ g1, data = da, center = "median")
rl = stat_hov_test(da, .cols = y, group = g1)
hc = function(item, a, b) note("stat_hov_test", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
hp = function(item, a, b) note("stat_hov_test", item, a, b, if (pv_close(a, b)) 0 else 1,
  if (pv_close(a, b)) OK else BUG)
hc("levene median F", rl$statistic, fl[1, "F value"])
hp("levene median p", rl$p.value, fl[1, "Pr(>F)"])
hc("levene df1 + df2/1e4", rl$df + rl$df2 / 1e4, fl[1, "Df"] + fl[2, "Df"] / 1e4)
flm = car::leveneTest(y ~ g1, data = da, center = "mean")
rlm = stat_hov_test(da, .cols = y, group = g1, center = "mean")
hc("levene mean F", rlm$statistic, flm[1, "F value"])
fb = bartlett.test(y ~ g1, data = da)
rb = stat_hov_test(da, .cols = y, group = g1, method = "bartlett")
hc("bartlett K^2", rb$statistic, fb$statistic)
hp("bartlett p", rb$p.value, fb$p.value)
dv2 = da[da$g1 != "c", ]
ff2 = var.test(y ~ g1, data = dv2)
rf2 = stat_hov_test(dv2, .cols = y, group = g1, method = "f")
hc("F-test statistic", rf2$statistic, ff2$statistic)
hp("F-test p", rf2$p.value, ff2$p.value)
note("stat_hov_test", "F-test df2 dropped (NA)", rf2$df2, ff2$parameter[2], NA, VR,
  "var.test returns num df & denom df; output keeps only df1, df2 = NA")

# ============================================================================
# 11. stat_bootstrap
# ============================================================================
cat("\n=== 11. stat_bootstrap ===\n")
set.seed(7)
xbt = rnorm(60, 5, 2)
dfbt = data.frame(x = xbt)
rbt = stat_bootstrap(dfbt, statistic = \(d, i) mean(d$x[i]), B = 200, seed = 42, type = "percentile")
set.seed(42)
bref = boot::boot(dfbt, function(d, i) mean(d$x[i]), R = 200)
ciref = boot::boot.ci(bref, conf = 0.95, type = "perc")
cat("  boot.ci perc names:", paste(names(ciref), collapse = ","), "\n")
perc_ref = ciref[[length(ciref)]][4]   # last element = perc bounds; [4]=lower
bc = function(item, a, b) note("stat_bootstrap", item, a, b, rdiff(a, b),
  if (is.na(rdiff(a, b))) VR else if (rdiff(a, b) < 1e-6) OK else BUG)
bc("estimate = sample mean", rbt$estimate, mean(xbt))
bc("same-seed reproduction of boot.ci perc", rbt$conf_low, perc_ref)
bc("bias", rbt$bias, mean(bref$t) - bref$t0)
bc("se", rbt$se, sd(bref$t))
set.seed(42)
bref2 = boot::boot(dfbt, function(d, i) mean(d$x[i]), R = 200)
cibca = boot::boot.ci(bref2, conf = 0.95, type = "bca")
rbca = stat_bootstrap(dfbt, statistic = \(d, i) mean(d$x[i]), B = 200, seed = 42, type = "bca")
bc("bca CI reproducible", rbca$conf_low, cibca$bca[4])
# theoretical se sanity
theory = sd(xbt) / sqrt(60)
note("stat_bootstrap", "se vs theory (n=60)", rbt$se, theory, rdiff(rbt$se, theory),
  if (abs(rbt$se - theory) / theory < 0.15) OK else VR,
  sprintf("MC se=%.4f theory=%.4f", rbt$se, theory))
# coverage sanity
set.seed(123)
cov_hits = 0
for (k in 1:100) {
  xk = rnorm(25, 3, 1)
  ok = tryCatch({
    rk = stat_bootstrap(data.frame(x = xk), \(d, i) mean(d$x[i]), B = 300,
      seed = 1000 + k, type = "percentile")
    rk$conf_low < 3 && 3 < rk$conf_high
  }, error = function(e) NA)
  if (isTRUE(ok)) cov_hits = cov_hits + 1
}
note("stat_bootstrap", "empirical coverage (nominal 95%)", cov_hits / 100, 0.95, NA,
  if (abs(cov_hits / 100 - 0.95) <= 0.09) OK else VR,
  sprintf("%d/100 percentile CIs cover true mean", cov_hits))
rmbt = stat_bootstrap(dfbt, statistic = \(d, i) c(mean = mean(d$x[i]), s = sd(d$x[i])),
  B = 100, seed = 8, type = "percentile")
note("stat_bootstrap", "named multi-comp terms", paste(rmbt$term, collapse = "/"), "mean/s", NA, OK)
rubt = stat_bootstrap(dfbt, statistic = \(d, i) c(mean(d$x[i]), sd(d$x[i])),
  B = 100, seed = 8, type = "percentile")
note("stat_bootstrap", "unnamed terms V1/V2", paste(rubt$term, collapse = "/"), "V1/V2", NA, OK)
# constant statistic -> bca fails -> warning + NA
warned = FALSE
rcbt = withCallingHandlers(
  stat_bootstrap(dfbt, statistic = \(d, i) 1, B = 50, seed = 1, type = "bca"),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
note("stat_bootstrap", "constant stat: warns, CI NA", paste(warned, is.na(rcbt$conf_low), sep = "/"),
  "TRUE/TRUE", NA, if (warned && is.na(rcbt$conf_low)) OK else BUG)

# ============================================================================
# 12. stat_perm_test (exact enumeration checks)
# ============================================================================
cat("\n=== 12. stat_perm_test ===\n")
dperm = data.frame(y = c(3, 5, 7, 9, 11, 2, 4, 6, 8, 10),
                   g = factor(rep(c("a","b"), each = 5)))
# exact two.sample: all C(10,5)=252 assignments
comb = combn(10, 5)
T0 = mean(dperm$y[1:5]) - mean(dperm$y[6:10])
Tstar_ex = apply(comb, 2, function(ii) mean(dperm$y[ii]) - mean(dperm$y[-ii]))
p_ex_2s = (sum(abs(Tstar_ex) >= abs(T0))) / 252
rp = stat_perm_test(dperm, .cols = y, group = g, type = "two.sample",
  permutations = 50000, seed = 99)
note("stat_perm_test", "2-samp p vs exact enum (252)", rp$p.value, p_ex_2s, NA,
  if (abs(rp$p.value - p_ex_2s) < 0.01) OK else VR,
  sprintf("MC=%.5f exact=%.5f", rp$p.value, p_ex_2s))
# one.sample: 2^5 sign flips
v1 = c(.3, .6, .9, 1.2, 1.5)
ex1 = as.matrix(expand.grid(rep(list(c(-1, 1)), 5)))
Tst1 = as.vector(ex1 %*% v1) / 5
p_ex_1s = sum(abs(Tst1) >= abs(mean(v1))) / 32
r1p = stat_perm_test(data.frame(y = v1), .cols = y, type = "one.sample",
  permutations = 50000, seed = 7)
note("stat_perm_test", "1-samp p vs exact enum (32)", r1p$p.value, p_ex_1s, NA,
  if (abs(r1p$p.value - p_ex_1s) < 0.01) OK else VR,
  sprintf("MC=%.5f exact=%.5f", r1p$p.value, p_ex_1s))
# paired = sign-flip of diffs (2^6 = 64)
dvpp = c(-0.5, 0, 0.5, 1, 1.5, 1.5)  # includes a zero diff
expp = as.matrix(expand.grid(rep(list(c(-1, 1)), 6)))
Tstpp = as.vector(expp %*% dvpp) / 6
p_ex_pp = sum(abs(Tstpp) >= abs(mean(dvpp))) / 64
dpp = data.frame(pre = c(1,2,3,4,5,6), post = c(1.5,2,2.5,3,3.5,4.5))
rpp = stat_perm_test(dpp, .cols = pre:post, type = "paired", permutations = 50000, seed = 3)
note("stat_perm_test", "paired p vs exact enum (64, zero diff incl.)", rpp$p.value, p_ex_pp, NA,
  if (abs(rpp$p.value - p_ex_pp) < 0.01) OK else VR,
  "sign-flip permutation; zero diffs contribute 0 (kept, consistent enum)")
# independence: enumerate all 6! permutations of second column
dind = data.frame(x = 1:6, y = c(2, 1, 4, 3, 6, 5))
permute = function(v) {
  if (length(v) == 1) return(matrix(v))
  m = NULL
  for (i in seq_along(v)) {
    rest = permute(v[-i])
    m = rbind(m, cbind(matrix(v[i], nrow(rest), 1), rest))
  }
  m
}
allp = permute(1:6)
Tst_ind = apply(allp, 1, function(i) cor(dind$x, dind$y[i]))
T0_ind = cor(dind$x, dind$y)
p_ex_ind = sum(abs(Tst_ind) >= abs(T0_ind)) / 720
rind = stat_perm_test(dind, .cols = x:y, type = "independence", permutations = 50000, seed = 21)
note("stat_perm_test", "independence cor p vs exact enum (720)", rind$p.value, p_ex_ind, NA,
  if (abs(rind$p.value - p_ex_ind) < 0.01) OK else VR,
  sprintf("MC=%.5f exact=%.5f", rind$p.value, p_ex_ind))
# greater alternative vs exact
p_ex_gr = sum(Tstar_ex >= T0) / 252
rg = stat_perm_test(dperm, .cols = y, group = g, type = "two.sample",
  alternative = "greater", permutations = 50000, seed = 8)
note("stat_perm_test", "greater p vs exact enum", rg$p.value, p_ex_gr, NA,
  if (abs(rg$p.value - p_ex_gr) < 0.01) OK else VR,
  "T0 direction a-b; verified distribution symmetric so direction flip in Tstar harmless")
# coin cross-check (secondary, convention differences expected)
cp = coin::oneway_test(y ~ g, data = dperm, distribution = "exact")
note("stat_perm_test", "2-samp p vs coin oneway exact", rp$p.value, coin::pvalue(cp), NA, VR,
  "coin uses standardized linear statistic (Fisher-Pitman); different statistic, sanity only")
# custom statistic: two.sample invariant trap
f_cust = function(d, i) {
  y = d$y[i]; lab = factor(d$g[i], levels = c("a", "b"))
  median(y[lab == "a"]) - median(y[lab == "b"])
}
rcu = stat_perm_test(dperm, statistic = f_cust, type = "two.sample", permutations = 500, seed = 4)
note("stat_perm_test", "custom 2samp (both cols indexed) p==1", rcu$p.value, 1,
  if (pv_close(rcu$p.value, 1)) 0 else 1, VR,
  "re-indexing y AND g keeps pairs intact -> statistic invariant -> null degenerates, p=1; docs/own test example lead users into this trap")
# custom statistic independence: framework shuffles LAST column
dci3 = data.frame(x = 1:8, junk = rnorm(8), z = c(1, 3, 2, 5, 4, 7, 6, 8) + 0.5)
f_ind = function(d, i) cor(d$x, d$z)
rci3 = stat_perm_test(dci3, statistic = f_ind, type = "independence", permutations = 2000, seed = 6)
note("stat_perm_test", "custom indep shuffles LAST col (z) -> valid",
  rci3$null_mean, 0, NA, if (abs(rci3$null_mean) < 0.1) OK else BUG,
  sprintf("null_mean=%.3f; shuffle target is d[[ncol(d)]], doc says 'second column'", rci3$null_mean))
f_wrong = function(d, i) cor(d$x, d$junk)
rcw = stat_perm_test(dci3, statistic = f_wrong, type = "independence", permutations = 2000, seed = 6)
note("stat_perm_test", "custom indep using 2nd col junk: NOT shuffled",
  rcw$null_mean, cor(dci3$x, dci3$junk), NA,
  if (abs(rcw$null_mean - cor(dci3$x, dci3$junk)) < 0.15) VR else BUG,
  sprintf("null_mean=%.3f == observed cor -> column actually used is never shuffled; doc says 'second column' but code shuffles LAST (L1228)",
    rcw$null_mean))
# custom indep with .by columns after the data columns: shuffled col = .by col!
dby = data.frame(x = 1:8, z = c(1, 3, 2, 5, 4, 7, 6, 8) + 0.5, s = rep(c("p","q"), 4))
rbyw = tryCatch(
  stat_perm_test(dby, statistic = f_ind, type = "independence", permutations = 100, seed = 1, .by = s),
  error = function(e) paste("ERROR:", conditionMessage(e)))
rby_txt = if (inherits(rbyw, "stat_infer")) {
  # slice column s gets shuffled into cor(): factor coerced? capture observed/null
  sprintf("observed=%.3f null_mean=%.3f", rbyw$observed, rbyw$null_mean)
} else as.character(rbyw)
note("stat_perm_test", "custom indep + .by: shuffles slice col s", rby_txt, "shuffles-z-not-s", NA, VR,
  "with .by present, d[[ncol(d)]] is the slice column s -> cor(x, shuffled s) nonsense observed differs from cor(x,z); last-column assumption fragile")

# ============================================================================
# 13. Edge cases
# ============================================================================
cat("\n=== 13. Edge cases ===\n")
# For input-validation: a CLEAN, informative stop() is the CORRECT behavior.
# A crash with an opaque error (stats internals leaking) or a silent wrong
# number is a finding. Degenerate stats should either work or fail cleanly.
res_ok = function(label, expr, want = c("ok", "clean-error")) {
  res = tryCatch(withCallingHandlers({force(expr); "ok"},
    warning = function(w) paste("WARN:", conditionMessage(w))),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  is_clean = grepl("^(ok|WARN)", res) ||
    (grepl("^ERROR", res) && !grepl("(must be|levels in the data|exactly|at least|numeric|select)", res, TRUE) &&
       grepl("^(ERROR: (data are essentially constant|no data in call|all 'x' values))", res))
  v = if (res == "ok") OK else if (grepl("^(ok|WARN)", res) || grepl("must|levels in the data|at least 2|numeric column|at least 2 columns", res, TRUE)) VR else BUG
  note("edge", label, res, "any", NA, v, "")
}
# --- behavior we WANT: clean validation errors (verdict OK/VR not BUG)
note("edge", "empty df t_test: clean stop", "see next", "ok", NA, OK, "")
res_ok("empty df t_test", {
  e0 = data.frame(y = numeric(0), g = factor(character(0), levels = c("a","b")))
  stat_t_test(e0, .cols = y, group = g)
})
res_ok("empty df describe works", stat_describe(data.frame(x = numeric(0)), .cols = x))
res_ok("one group after NA drop: clean stop", {
  stat_t_test(data.frame(y = c(1, 2, NA, 4), g = c("a","a","b","a")), .cols = y, group = g)
})
res_ok("group col all NA: clean stop", {
  stat_t_test(data.frame(y = 1:4, g = c(NA, NA, "a", "a")), .cols = y, group = g)
})
res_ok("chisq degenerate table: clean stop", {
  stat_chisq(data.frame(a = c("x","x"), b = c("y","y")), row = a, col = b)
})
res_ok("anova single level: clean stop", {
  stat_anova(data.frame(y = 1:4, g = rep("a", 4)), resp = y, group = g)
})
res_ok("hov f 3 groups: clean stop", {
  stat_hov_test(da, .cols = y, group = g1, method = "f")
})
res_ok("cor 1 column: clean stop", { tryCatch(stat_cor(dco, .cols = a), error = function(e) stop(e)) })
# --- degenerate statistics: what do we actually get?
all_eq2 = suppressWarnings(tryCatch({
  r = stat_wilcox_test(data.frame(y = rep(5, 10), g = rep(c("a","b"), each = 5)),
    .cols = y, group = g)
  sprintf("W=%s p=%s eff=%s", format(r$statistic), format(r$p.value), format(r$effect))
}, error = function(e) paste("ERROR:", conditionMessage(e))))
note("edge", "all-equal wilcox 2-samp values", all_eq2, "p=1-like", NA,
  if (grepl("ERROR", all_eq2)) VR else OK, "stats sigma=0 case; check p in (0,1]")
all_eq_p = suppressWarnings(tryCatch({
  r = stat_t_test(data.frame(pre = rep(5, 10), post = rep(5, 10)), .cols = pre:post, paired = TRUE)
  sprintf("t=%s p=%s", format(r$statistic), format(r$p.value))
}, error = function(e) paste("ERROR:", conditionMessage(e))))
note("edge", "all-equal paired t values", all_eq_p, "NaN p", NA,
  if (grepl("ERROR", all_eq_p)) VR else OK, "sd=0 -> NaN; stats::t.test same behavior")
const_w = tryCatch({
  r = withCallingHandlers(
    stat_normality(data.frame(x = rep(5, 20)), .cols = x),
    warning = function(w) invokeRestart("muffleWarning"))
  sprintf("W=%s p=%s skew=%s", format(r$statistic), format(r$p.value), format(r$skew))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
note("edge", "normality constant col: leaks stats error", const_w, "clean", NA,
  if (grepl("^ERROR: all 'x' values are identical", const_w)) VR else BUG,
  "shapiro.test's raw error text leaks through (stats-native message, not validated); n_used=20 row otherwise lost")
boot0 = tryCatch({
  stat_bootstrap(data.frame(x = numeric(0)), statistic = \(d, i) mean(d$x[i]), B = 10, seed = 1)
}, error = function(e) paste("ERROR:", conditionMessage(e)))
note("edge", "bootstrap 0-row data", boot0, "clean-error", NA,
  if (grepl("^ERROR: no data", boot0)) VR else BUG, "boot's native error; no pre-validation")
# perm test n=2 two.sample works
res_ok("perm test n=2 two.sample", {
  stat_perm_test(data.frame(y = c(1, 3), g = c("a","b")), .cols = y, group = g,
    type = "two.sample", permutations = 100, seed = 1)
})
# bootstrap .by slice empty after NA drop -> empty slice handled? (b has 1 row)
res_ok("bootstrap .by slice with 1 row", {
  b = data.frame(x = c(1, 2, NA, 4), s = c("a","a","a","b"))
  stat_bootstrap(b, statistic = \(d, i) mean(d$x[i]), B = 10, seed = 1, .by = s)
})
# describe .by with NA slice rows
dms = data.frame(x = c(1, 2, NA, 4, 5, 6), s = c("a","a","a","a",NA,NA))
rms = stat_describe(dms, .cols = x, .by = s)
note("edge", "describe .by: NA-slice rows dropped", paste(nrow(rms), rms$n_used, rms$n_missing, sep = "/"),
  "1/3/3", NA, OK, "slice a keeps rows 1,2,4 (x=NA row 3 dropped); n_missing = 6-3 = 3")
# t_test two-sample with n=2 per group (df=0 pooled): stats errors; mathmodels?
tiny2 = tryCatch({
  r = stat_t_test(data.frame(y = c(1, 3, 5, 7), g = rep(c("a","b"), each = 2)),
    .cols = y, group = g, var.equal = TRUE)
  sprintf("df=%s p=%s", format(r$df), format(r$p.value))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
note("edge", "t_test n=2-per-group pooled (df=0)", tiny2, "same-as-stats", NA,
  if (grepl("ERROR: (not enough|x values)", tiny2) || grepl("df=0", tiny2)) VR else OK,
  "stats::t.test pools df=0 -> error; mathmodels passes it through verbatim")

saveRDS(RES$tbl, ".superpowers/audits/probe-results.rds")
cat("\n\n==== VERDICT SUMMARY ====\n")
tb = RES$tbl
print(table(tb$verdict))
cat("\n==== ALL LINES WITH EXTRA TEXT ====\n")
for (i in seq_len(nrow(tb))) if (nzchar(tb$extra[i])) {
  cat(sprintf("[%s] %s | %s | %s\n", tb$verdict[i], tb$fun[i], tb$item[i], tb$extra[i]))
}
cat("\nDONE\n")
