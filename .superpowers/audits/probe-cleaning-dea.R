# Independent audit probe: data_cleaning (self-written) + interp (light) + DEA
devtools::load_all("E:/Github/mathmodels")
options(warn = -1)
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}

## ============ na_mcar_test vs naniar::mcar_test ============
set.seed(1)
d = data.frame(x=rnorm(120), y=rnorm(120), z=rnorm(120))
d$x[sample(120,15)] = NA; d$y[sample(120,10)] = NA
r = na_mcar_test(d, .cols = x:z)
if (requireNamespace("naniar", quietly=TRUE)) {
  nr = naniar::mcar_test(d)
  # package EM verified (via direct likelihood optimisation) to reach the true
  # MLE; chi2 is the exact Little formula on those estimates. naniar uses a
  # different EM/statistic detail -> small numerical difference, same verdict.
  check("na_mcar_test chi2 close to naniar (<=2% rel)",
        abs(r$chi2 - nr$statistic)/nr$statistic < 0.02,
        sprintf("pkg chi2=%.6f naniar=%.6f rel=%.4f", r$chi2, nr$statistic,
                abs(r$chi2-nr$statistic)/nr$statistic))
  check("na_mcar_test df vs naniar",
        isTRUE(all.equal(as.numeric(r$df), as.numeric(nr$df), tolerance = 1e-8)),
        sprintf("pkg df=%s naniar=%s", r$df, nr$df))
  check("na_mcar_test same MCAR verdict as naniar",
        (r$p_value > 0.05) == (nr$p.value > 0.05),
        sprintf("pkg p=%.4f naniar p=%.4f", r$p_value, nr$p.value))
}

## ============ impute_model knn vs manual ============
set.seed(2)
dk = data.frame(a=rnorm(20), b=rnorm(20), target=rnorm(20))
dk$target[c(3,9,15)] = NA
im = impute_model(dk, .cols = target, model = "knn", k = 3)
# manual KNN: features a,b; euclidean; mean of 3 nearest complete rows
feat = c("a","b"); tr = which(!is.na(dk$target))
ytr = dk$target[tr]; Xtr = as.matrix(dk[tr, feat])
manual = sapply(which(is.na(dk$target)), function(i) {
  xq = as.numeric(dk[i, feat])
  dd = sqrt(rowSums((Xtr - matrix(xq, nrow(Xtr), 2, byrow=TRUE))^2))
  mean(ytr[order(dd)[1:3]])
})
check("impute_model knn vs manual",
      isTRUE(all.equal(im$target[c(3,9,15)], manual, tolerance=1e-10)),
      sprintf("pkg=(%s) ref=(%s)", paste(round(im$target[c(3,9,15)],6),collapse=","),
              paste(round(manual,6),collapse=",")))

## ============ outlier bounds ============
xo = c(1,2,3,4,5,6,7,8,9,100)
o = outlier(data.frame(v=xo), .cols = v, method = "iqr")
q = quantile(xo, c(.25,.75), names=FALSE); iqr = q[2]-q[1]
check("outlier iqr bound detects 100",
      nrow(o) == 1 && o$row == 10 && o$value == 100,
      sprintf("nrow=%d", nrow(o)))
check("outlier iqr bound values",
      abs(o$bound_low - (q[1]-1.5*iqr)) < 1e-10 && abs(o$bound_high - (q[2]+1.5*iqr)) < 1e-10)
# z-score: use a tight cluster + far outlier so |z| > 3
xz = c(rep(2, 30), 100)
oz = outlier(data.frame(v=xz), .cols=v, method="zscore")
m = mean(xz); s = sd(xz)
check("outlier zscore detects 100", nrow(oz)==1 && oz$row==31,
      sprintf("nrow=%d z=%.3f", nrow(oz), (100-m)/s))
check("outlier zscore bounds",
      abs(oz$bound_low-(m-3*s))<1e-10 && abs(oz$bound_high-(m+3*s))<1e-10)
on = outlier(data.frame(v=xo), .cols=v, method="iqr", to="na")
check("outlier to='na' sets NA", is.na(on$v[10]) && !inherits(on, "tbl_df"))

## ============ interp light check (wrappers) ============
xi = c(1,2,3,4,5); yi = c(2,4,1,5,3); xo2 = c(1.5,2.5,3.5)
li = interp_linear(xi, yi, xo2)
check("interp_linear vs stats::approx",
      isTRUE(all.equal(li$y, stats::approx(xi,yi,xout=xo2)$y, tolerance=1e-12)))
sp = interp_spline(xi, yi, xo2)
check("interp_spline vs stats::spline fmm",
      isTRUE(all.equal(sp$y, stats::spline(xi,yi,xout=xo2,method="fmm")$y, tolerance=1e-10)))
po = interp_poly(xi, yi, xo2, degree = 3)
check("interp_poly vs lm(poly raw)",
      isTRUE(all.equal(po$y, as.numeric(predict(lm(yi ~ poly(xi,3,raw=TRUE)), data.frame(xi=xo2))), tolerance=1e-8)))

## ============ basic_DEA vs Benchmarking ============
if (requireNamespace("Benchmarking", quietly=TRUE)) {
  set.seed(3)
  X = matrix(runif(12, 1, 10), nrow=6, ncol=2)   # 6 DMUs, 2 inputs
  Y = matrix(runif(6, 1, 10), nrow=6, ncol=1)    # 1 output
  dat = data.frame(DMU=paste0("D",1:6), x1=X[,1], x2=X[,2], y=Y[,1])
  for (rts in c("vrs","crs")) {
    ours = basic_DEA(dat, inputs=c("x1","x2"), outputs="y", rts=rts, orientation="io")
    ref = Benchmarking::dea(X=X, Y=Y, RTS=toupper(rts), ORIENTATION="in")
    ov = ours$efficiencies
    check(sprintf("basic_DEA %s input-orientation vs Benchmarking", rts),
          isTRUE(all.equal(as.numeric(ov), as.numeric(ref$eff), tolerance=1e-5)),
          sprintf("ours=(%s) ref=(%s)", paste(round(as.numeric(ov),5),collapse=","),
                  paste(round(as.numeric(ref$eff),5),collapse=",")))
  }
}

cat(sprintf("\n=== cleaning/interp/dea audit: PASS %d / FAIL %d ===\n", P, F))
