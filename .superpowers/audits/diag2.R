devtools::load_all("E:/Github/mathmodels")
options(warn = -1)

## DEA: full return names + efficiency extraction
set.seed(3)
X = matrix(runif(12, 1, 10), nrow=6, ncol=2)
Y = matrix(runif(6, 1, 10), nrow=6, ncol=1)
dat = data.frame(DMU=paste0("D",1:6), x1=X[,1], x2=X[,2], y=Y[,1])
o = basic_DEA(dat, inputs=c("x1","x2"), outputs="y", rts="crs", orientation="io")
cat("NAMES:", paste(names(o), collapse=", "), "\n")
str(o, max.level=1)

## na_mcar: self-consistency of the Little statistic with the package EM
set.seed(1)
d = data.frame(x=rnorm(120), y=rnorm(120), z=rnorm(120))
d$x[sample(120,15)] = NA; d$y[sample(120,10)] = NA
Xm = as.matrix(d); em = mathmodels:::.em_norm_mle(Xm)
mu = em$mu; Sig = em$sigma; Si = solve(Sig)
miss = is.na(Xm); key = apply(miss, 1, paste, collapse="")
grp = split(seq_len(nrow(Xm)), key)
chi2_manual = 0; df_manual = 0
for (g in grp) {
  mp = miss[g[1],]; oc = which(!mp)
  mg = colMeans(Xm[g, oc, drop=FALSE]); dif = mg - mu[oc]
  chi2_manual = chi2_manual + length(g) * as.numeric(t(dif) %*% Si[oc,oc,drop=FALSE] %*% dif)
  df_manual = df_manual + length(oc)
}
df_manual = df_manual - ncol(Xm)
cat(sprintf("manual chi2=%.8f df=%d vs pkg %.8f df=%d\n", chi2_manual, df_manual,
            na_mcar_test(d, .cols=x:z)$chi2, na_mcar_test(d, .cols=x:z)$df))
