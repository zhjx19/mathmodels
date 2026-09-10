# Independent audit probe: pred_grey.R (GM11, GM1N, DGM21, verhulst, combine_preds)
# Reference implementations coded from textbook formulas, NOT package code.

devtools::load_all("E:/Github/mathmodels")
set.seed(2026)
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (ok) { P = P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F = F + 1; cat("FAIL:", name, "-", detail, "\n") }
  assign("P", P, envir = globalenv()); assign("F", F, envir = globalenv())
}

## ---------------------------------------------------------------------------
## 1. GM(1,1)
##    Textbook: X1 = cumsum(X0); z(k) = 0.5*(X1(k)+X1(k+1)); Y = X0(2..n)
##    [a,b] = (B'B)^-1 B'Y, B = [-z, 1]
##    x1_hat(k) = (X0(1) - b/a) e^{-a(k-1)} + b/a     (k = 1..n+1)
##    x0_hat(k) = x1_hat(k) - x1_hat(k-1)
## ---------------------------------------------------------------------------
gm11_ref = function(X) {
  n = length(X)
  X1 = cumsum(X)
  z = 0.5 * (X1[-n] + X1[-1])
  B = cbind(-z, 1)
  Y = X[-1]
  u = solve(t(B) %*% B, t(B) %*% Y)
  a = u[1]; b = u[2]
  X1h = (X[1] - b/a) * exp(-a * (0:n)) + b/a   # k = 1..n+1
  x0h = c(X1h[1], diff(X1h))
  list(a = a, b = b, fitted = x0h[1:n], pnext = x0h[n+1])
}

x_ser = c(100, 120, 145, 175, 210)   # doc example
r  = GM11(x_ser)
ref = gm11_ref(x_ser)
check("GM11 param a", isTRUE(all.equal(r$u[1], ref$a, tolerance = 1e-8)),
      sprintf("pkg=%.8f ref=%.8f", r$u[1], ref$a))
check("GM11 param b", isTRUE(all.equal(r$u[2], ref$b, tolerance = 1e-8)),
      sprintf("pkg=%.8f ref=%.8f", r$u[2], ref$b))
check("GM11 fitted", isTRUE(all.equal(as.numeric(r$fitted), ref$fitted, tolerance = 1e-6)),
      sprintf("maxdiff=%.3e", max(abs(as.numeric(r$fitted) - ref$fitted))))
check("GM11 pnext", isTRUE(all.equal(as.numeric(r$pnext), ref$pnext, tolerance = 1e-6)),
      sprintf("pkg=%.4f ref=%.4f", r$pnext, ref$pnext))
check("GM11 f(n+1) == pnext", isTRUE(all.equal(r$f(6), r$pnext, tolerance = 1e-8)))

# exact recovery on pure exponential
Xe = 100 * exp(0.05 * (0:7))
re = GM11(Xe)
check("GM11 exact recovery on exp data",
      max(abs(re$fitted - Xe)) / 100 < 1e-3,
      sprintf("maxrelerr=%.2e", max(abs(re$fitted - Xe)) / 100))

# edges
check("GM11 constant series stops", inherits(tryCatch(GM11(rep(5,6)), error = \(e) e), "error"))
check("GM11 n=3 stops", inherits(tryCatch(GM11(c(1,2,3)), error = \(e) e), "error"))
check("GM11 negative stops", inherits(tryCatch(GM11(c(-1,2,3,4,5)), error = \(e) e), "error"))

## ---------------------------------------------------------------------------
## 2. DGM(2,1): package regresses diff(X) on [-X(2..n), 1] via OLS:
##    X(k+1) - X(k) = -a X(k) + b  =>  X(k+1) = (1-a) X(k) + b   [AR(1)]
##    Verify params + that fitted/pnext follow the AR(1) closed form
##    X_hat(k) = (x1 - b/a)(1-a)^{k-1} + b/a
## ---------------------------------------------------------------------------
dgm21_ref = function(X) {
  n = length(X)
  # Liu Sifeng DGM(2,1): [x0(k)-x0(k-1)] + a x0(k) = b, k=2..n
  # => diff(X)[i] paired with X[i+1]
  B = cbind(-X[-1], 1)
  Y = diff(X)
  u = solve(t(B) %*% B, t(B) %*% Y)
  a = u[1]; b = u[2]
  x1_0 = X[1]
  # time response of x1 (AGO): x1(k) = (b/a^2 - x0/a)e^{-a(k-1)}
  #   + (b/a)(k-1) + (1+a)/a x0 - b/a^2 ; fitted = IAGO
  fk = function(k) (b/a^2 - x1_0/a)*exp(-a*(k-1)) + (b/a)*(k-1) +
    (1+a)/a*x1_0 - b/a^2
  X1h = fk(1:(n+1))
  x0h = c(X1h[1], diff(X1h))
  list(a = a, b = b, fitted = x0h[1:n], pnext = x0h[n+1])
}
x_d = c(2.874, 3.278, 3.39, 3.679, 3.77, 3.8)
rd = DGM21(x_d)
rd_ref = dgm21_ref(x_d)
check("DGM21 params match independent OLS",
      isTRUE(all.equal(as.numeric(rd$u), c(rd_ref$a, rd_ref$b), tolerance = 1e-8)),
      sprintf("pkg a=%.6f b=%.6f; ref a=%.6f b=%.6f", rd$u[1], rd$u[2], rd_ref$a, rd_ref$b))
# package f(k) is claimed to be the X1 (cumulated) response; IAGO via diff.
# Check consistency: fitted == IAGO of f, and pnext == f(n+1)-f(n)
fh = rd$f(1:(length(x_d)+1))
fit2 = c(fh[1], diff(fh))
check("DGM21 fitted = IAGO(f)", isTRUE(all.equal(as.numeric(rd$fitted), fit2[1:6], tolerance = 1e-6)),
      sprintf("maxdiff=%.3e", max(abs(as.numeric(rd$fitted) - fit2[1:6]))))
check("DGM21 pnext = f(n+1)-f(n)",
      isTRUE(all.equal(as.numeric(rd$pnext), fit2[7], tolerance = 1e-6)))
# closed-form: package fitted must match independent reconstruction of the
# textbook DGM(2,1) time response (params from independent OLS)
check("DGM21 fitted = independent closed-form reconstruction",
      max(abs(as.numeric(rd$fitted) - rd_ref$fitted)) < 1e-6,
      sprintf("maxdiff=%.3e", max(abs(as.numeric(rd$fitted) - rd_ref$fitted))))
check("DGM21 pnext = independent closed form",
      isTRUE(all.equal(as.numeric(rd$pnext), rd_ref$pnext, tolerance = 1e-6)))
# discretization-bias note (grey-model inherent, NOT a bug): on data from the
# DGM closed form with (a,b), the boundary-differenced pair biases the OLS;
# dropping it recovers the theoretical discrete relation.
a_t = 0.2; b_t = 1.5; x1_0 = 3
fk_t = function(k) (b_t/a_t^2 - x1_0/a_t)*exp(-a_t*(k-1)) + (b_t/a_t)*(k-1) +
  (1+a_t)/a_t*x1_0 - b_t/a_t^2
X1_t = fk_t(1:9); xr = c(X1_t[1], diff(X1_t))
Bx = cbind(-xr[-1], 1); Yx = diff(xr)
u6 = solve(t(Bx[-1, , drop = FALSE]) %*% Bx[-1, , drop = FALSE],
           t(Bx[-1, , drop = FALSE]) %*% Yx[-1])
check("DGM21 discrete relation theory (drop boundary pair)",
      isTRUE(all.equal(as.numeric(u6), c(exp(a_t) - 1, b_t/a_t*(exp(a_t) - 1)),
                       tolerance = 1e-6)),
      sprintf("recovered a=%.5f (theory e^a-1=%.5f)", u6[1], exp(a_t) - 1))
# fit quality on the documented example (modest relative error expected)
check("DGM21 fit quality on doc example",
      mean(abs(as.numeric(rd$fitted) - x_d) / x_d) < 0.05,
      sprintf("mean relerr=%.3e", mean(abs(as.numeric(rd$fitted) - x_d) / x_d)))

## ---------------------------------------------------------------------------
## 3. verhulst: on X1 = cumsum(X): dX1/dt + a X1 = b X1^2
##    OLS: Y = X0(2..n) on B = [-Z, Z^2] (Z = mean-adjacent X1)
##    X1_hat(k) = a*X1(1) / (b*X1(1) + (a - b*X1(1)) e^{a(k-1)})
## ---------------------------------------------------------------------------
verhulst_ref = function(X) {
  n = length(X)
  X1 = cumsum(X)
  Z = 0.5 * (X1[-n] + X1[-1])
  B = cbind(-Z, Z^2)
  u = solve(t(B) %*% B, t(B) %*% X[-1])
  a = u[1]; b = u[2]
  X1h = a * X[1] / (b * X[1] + (a - b * X[1]) * exp(a * (0:n)))
  x0h = c(X1h[1], diff(X1h))
  list(a = a, b = b, fitted = x0h[1:n], pnext = x0h[n+1])
}
x_v = c(4.93, 2.33, 3.87, 4.35, 6.63, 7.15, 5.37, 6.39, 7.81, 8.35)
rv = verhulst(x_v)
rv_ref = verhulst_ref(x_v)
check("verhulst params match independent OLS",
      isTRUE(all.equal(as.numeric(rv$u), c(rv_ref$a, rv_ref$b), tolerance = 1e-7)),
      sprintf("pkg a=%.6f b=%.3e; ref a=%.6f b=%.3e", rv$u[1], rv$u[2], rv_ref$a, rv_ref$b))
check("verhulst fitted match ref closed form",
      isTRUE(all.equal(as.numeric(rv$fitted), rv_ref$fitted, tolerance = 1e-5)),
      sprintf("maxdiff=%.3e", max(abs(as.numeric(rv$fitted) - rv_ref$fitted))))
check("verhulst pnext match ref",
      isTRUE(all.equal(as.numeric(rv$pnext), rv_ref$pnext, tolerance = 1e-4)))
# fit quality: Verhulst is an S-curve model; the doc series is non-monotone,
# so only a loose sanity bound is meaningful (implementation already verified
# exact against independent OLS + closed form).
check("verhulst fit quality sanity on doc example",
      mean(abs(as.numeric(rv$fitted) - x_v) / x_v) < 0.25,
      sprintf("mean relerr=%.3e", mean(abs(as.numeric(rv$fitted) - x_v) / x_v)))
# recovers near-true params on data from the Verhulst closed form (small bias)
a_v = 0.5; b_v = 0.02; X1_0 = 5
fv = function(k) a_v * X1_0 / (b_v * X1_0 + (a_v - b_v * X1_0) * exp(a_v * (k - 1)))
X1_v = fv(1:12); x0l = c(X1_v[1], diff(X1_v))
rv2 = verhulst(x0l)
check("verhulst near-true param recovery on closed-form data",
      abs(rv2$u[1] - a_v) < 0.02 && abs(rv2$u[2] - b_v) < 0.002,
      sprintf("a=%.5f (true %.2f) b=%.6f (true %.3f)", rv2$u[1], a_v, rv2$u[2], b_v))

## ---------------------------------------------------------------------------
## 4. GM(1,N): dX1_c/dt + a X1_c = sum b_i X1_i
##    OLS: Y = X0_c(2..n) on [-Z_c, X1_factors(2..n)]
##    X1c_hat(k) = (X1c(1) - (1/a) sum b_i X1_i(k)) e^{-a(k-1)} + (1/a) sum b_i X1_i(k)
## ---------------------------------------------------------------------------
gm1n_ref = function(dat) {
  X0 = as.matrix(dat); n = nrow(X0); m = ncol(X0)
  X1 = apply(X0, 2, cumsum)
  Z = -0.5 * (X1[-n, m] + X1[-1, m])
  B = cbind(Z, X1[-1, -m])
  u = solve(t(B) %*% B, t(B) %*% X0[-1, m])
  a = u[1]; bb = u[-1]
  X1h = sapply(1:n, \(k)
    (X1[1, m] - (1/a) * sum(X1[k, -m] * bb)) * exp(-a * (k - 1)) +
    (1/a) * sum(X1[k, -m] * bb))
  fitted = c(X1h[1], diff(X1h))
  list(fitted = fitted)
}
dat = data.frame(
  factor1 = c(50, 55, 60, 65, 70),
  factor2 = c(20, 22, 25, 28, 30),
  output  = c(100, 120, 145, 175, 210)
)
rn = GM1N(dat)
rn_ref = gm1n_ref(dat)
cat("GM1N fitted attrs: pkg names=", !is.null(names(rn$fitted)),
    " ref names=", !is.null(names(rn_ref$fitted)),
    " pkg class=", class(rn$fitted)[1], "\n", sep = "")
check("GM1N fitted match independent closed form",
      max(abs(as.numeric(rn$fitted) - as.numeric(rn_ref$fitted))) < 1e-6,
      sprintf("maxdiff=%.3e", max(abs(as.numeric(rn$fitted) - as.numeric(rn_ref$fitted)))))
S1 = sd(dat$output); resid = dat$output - rn$fitted
C_ref = sd(resid) / S1
P_ref = mean(abs(resid - mean(resid)) < 0.6745 * S1)
check("GM1N C metric", isTRUE(all.equal(rn$C, C_ref, tolerance = 1e-10)))
check("GM1N P metric", isTRUE(all.equal(rn$P, P_ref)))
pn = GM1N(dat, new_data = data.frame(factor1 = 75, factor2 = 32))$pred
check("GM1N new_data -> non-null pred", !is.null(pn), sprintf("pred=%.3f", pn))

## ---------------------------------------------------------------------------
## 5. combine_preds
## ---------------------------------------------------------------------------
cp = combine_preds(c(100, 102, 98))
x5 = c(100, 102, 98)
d5 = c(abs(100-102), abs(100-98), abs(102-98)); maxd = max(d5)
idx = expand.grid(1:3, 1:3)
R = matrix(cos(pi * (x5[idx[,1]] - x5[idx[,2]]) / (2 * maxd)), nrow = 3, byrow = TRUE)
w_ref = eigen(R)$vectors[, 1]; w_ref = w_ref / sum(w_ref)
check("combine_preds w = principal eigenvector (normalized)",
      isTRUE(all.equal(cp$w, w_ref, tolerance = 1e-10)),
      sprintf("pkg=(%s)", paste(round(cp$w, 6), collapse = ",")))
check("combine_preds a = w-weighted mean",
      isTRUE(all.equal(cp$a, sum(x5 * cp$w), tolerance = 1e-10)))
c2 = combine_preds(c(1, 3))
check("combine_preds n=2 equal weights",
      isTRUE(all.equal(c2$w, c(0.5, 0.5))) && isTRUE(all.equal(c2$a, 2)))
c3 = combine_preds(rep(7, 4))
check("combine_preds identical -> equal weights",
      isTRUE(all.equal(c3$w, rep(0.25, 4))) && isTRUE(all.equal(c3$a, 7)))
cp5 = combine_preds(c(10, 10.0000001, 20))
check("combine_preds weights non-negative",
      all(cp5$w >= -1e-12), sprintf("w=(%s)", paste(round(cp5$w, 4), collapse = ",")))
# negative eigenvalue sign flip safety: bimodal far apart
cp6 = combine_preds(c(0, 100, 100.0001))
check("combine_preds bimodal weights sane (non-negative)",
      all(cp6$w >= -1e-12), sprintf("w=(%s)", paste(round(cp6$w, 4), collapse = ",")))

cat(sprintf("\n=== pred_grey audit: PASS %d / FAIL %d ===\n", P, F))
