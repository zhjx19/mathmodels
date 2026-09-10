devtools::load_all("E:/Github/mathmodels")
X = data.frame(x1=c(10,20,15), x2=c(5,10,8))
cat("weight_cv names?", !is.null(names(weight_cv(X))), "\n")
Xe = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
cat("weight_entropy names?", !is.null(names(weight_entropy(Xe, c("+","-"))$w)), "\n")
Xc = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
cat("weight_critic names?", !is.null(names(weight_critic(Xc, c("+","-"))$w)), "\n")
set.seed(1); Xp = data.frame(a=rnorm(30), b=rnorm(30), c=rnorm(30), d=rnorm(30))
rp = weight_pca(Xp, nfs = 2)
pc_nostd = prcomp(Xp, center=FALSE, scale.=FALSE)
cat("pca lambda == no-std prcomp sdev^2 [1:2]?",
    isTRUE(all.equal(unname(rp$lambda), unname((pc_nostd$sdev^2)[1:2]), tolerance=1e-8)), "\n")
# with proper index (all "+") -> standardized then prcomp
rp2 = weight_pca(Xp, index = rep("+", 4), nfs = 2)
pc_std = prcomp(scale(Xp), center=FALSE, scale.=FALSE)
cat("pca lambda == std prcomp sdev^2 [1:2] (index='+')?",
    isTRUE(all.equal(unname(rp2$lambda), unname((pc_std$sdev^2)[1:2]), tolerance=1e-8)), "\n")
