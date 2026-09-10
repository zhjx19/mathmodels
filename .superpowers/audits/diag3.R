devtools::load_all("E:/Github/mathmodels")
set.seed(1)
d = data.frame(x=rnorm(120), y=rnorm(120), z=rnorm(120))
d$x[sample(120,15)] = NA; d$y[sample(120,10)] = NA
X = as.matrix(d)
em = mathmodels:::.em_norm_mle(X)

# observed-data log-likelihood
loglik = function(par) {
  mu = par[1:3]
  L = matrix(0, 3, 3); L[lower.tri(L, diag=TRUE)] = par[4:9]
  Sig = L %*% t(L)
  ll = 0
  for (i in seq_len(nrow(X))) {
    oc = which(!is.na(X[i,]))
    if (length(oc) == 0) next
    ll = ll + mvtnorm_dens(X[i,oc], mu[oc], Sig[oc,oc,drop=FALSE])
  }
  -ll
}
mvtnorm_dens = function(x, mu, S) {
  k = length(x); d = x - mu
  -0.5*(k*log(2*pi) + as.numeric(determinant(S, logarithm=TRUE)$modulus) +
        as.numeric(t(d) %*% solve(S, d)))
}
# starting values
p0 = c(em$mu, em$sigma[lower.tri(em$sigma, diag=TRUE)])
opt = optim(p0, loglik, method="BFGS", control=list(maxit=5000, reltol=1e-12))
mu_opt = opt$par[1:3]
L = matrix(0,3,3); L[lower.tri(L,diag=TRUE)] = opt$par[4:9]; Sig_opt = L %*% t(L)
cat("EM mu   :", round(em$mu, 6), "\n")
cat("optim mu:", round(mu_opt, 6), "\n")
cat("max |mu_EM - mu_optim| =", max(abs(em$mu - mu_opt)), "\n")
cat("max |Sig_EM - Sig_optim| =", max(abs(em$sigma - Sig_opt)), "\n")
cat("loglik EM =", -loglik(p0), " | loglik optim =", -opt$value, "\n")
