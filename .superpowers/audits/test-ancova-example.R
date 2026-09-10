devtools::load_all("E:/Github/mathmodels")
set.seed(1)
d = data.frame(
  g = rep(c("a", "b"), each = 10),
  x = rep(1:10, 2),
  y = 2 * rep(1:10, 2) + rep(c(0, 5), each = 10) + rnorm(20)
)
print(stat_ancova(d, resp = y, group = g, covariates = x))
cat("OK\n")
