devtools::load_all("E:/Github/mathmodels")
options(warn = -1)

## DEA return structure
set.seed(3)
X = matrix(runif(12, 1, 10), nrow=6, ncol=2)
Y = matrix(runif(6, 1, 10), nrow=6, ncol=1)
dat = data.frame(DMU=paste0("D",1:6), x1=X[,1], x2=X[,2], y=Y[,1])
o = basic_DEA(dat, inputs=c("x1","x2"), outputs="y", rts="crs", orientation="io")
cat("class:", class(o), "\n"); print(names(o)); print(o)
ref = Benchmarking::dea(X=X, Y=Y, RTS="CRS", ORIENTATION="in")
cat("benchmarking eff:", ref$eff, "\n")

## na_mcar: is the chi2 gap a convergence/estimator artifact?
set.seed(1)
d = data.frame(x=rnorm(120), y=rnorm(120), z=rnorm(120))
d$x[sample(120,15)] = NA; d$y[sample(120,10)] = NA
r1 = mathmodels:::na_mcar_test(d, .cols=x:z, maxit=200, tol=1e-8)$chi2
r2 = mathmodels:::na_mcar_test(d, .cols=x:z, maxit=5000, tol=1e-12)$chi2
cat(sprintf("chi2 tol1e-8: %.8f | tol1e-12: %.8f | naniar: %.8f\n", r1, r2, naniar::mcar_test(d)$statistic))
# compare EM mu/sigma conventions on a complete matrix
Xc = as.matrix(d)
# package EM
em = mathmodels:::.em_norm_mle(Xc)
cat("EM mu:", round(em$mu,6), "\n")
cat("colMeans(complete):", round(colMeans(Xc, na.rm=TRUE),6), "\n")
# naniar's internal: check if it uses cov * (n-1)/n or cov
cat("EM sigma diag:", round(diag(em$sigma),6), "\n")
cat("cov(n-1) diag:", round(diag(stats::cov(Xc, use="pairwise.complete.obs")),6), "\n")
