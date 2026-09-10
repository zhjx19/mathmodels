devtools::load_all("E:/Github/mathmodels")
set.seed(1)
d = data.frame(x = rnorm(40), z = rnorm(40))
d$y = 0.5 * d$x + rnorm(40)
d$x[c(3, 9)] = NA
imps = impute_multiple(d, .cols = x, m = 5, maxit = 3, seed = 2)
fits = lapply(imps, \(dat) lm(y ~ x, dat))
pool_fit(fits)
