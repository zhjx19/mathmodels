# Independent audit probe: eval_weights.R
devtools::load_all("E:/Github/mathmodels")
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}

## ---------------------------------------------------------------- AHP
A = matrix(c(1,1/2,4,3,3, 2,1,7,5,5, 1/4,1/7,1,1/2,1/3, 1/3,1/5,2,1,1, 1/3,1/5,3,1,1),
           byrow = TRUE, nrow = 5)
r = AHP(A)
# independent: principal eigenvector of A
ev = eigen(A); Lmax_ref = Re(ev$values[which.max(Mod(ev$values))])
w_ref = Re(ev$vectors[, which.max(Mod(ev$values))]); w_ref = w_ref / sum(w_ref)
check("AHP w = principal eigenvector",
      isTRUE(all.equal(r$w, w_ref, tolerance = 1e-10)),
      sprintf("pkg=(%s)", paste(round(r$w, 5), collapse = ",")))
check("AHP Lmax", isTRUE(all.equal(r$Lmax, Lmax_ref, tolerance = 1e-10)),
      sprintf("pkg=%.8f ref=%.8f", r$Lmax, Lmax_ref))
check("AHP CI = (Lmax-n)/(n-1)",
      isTRUE(all.equal(r$CI, (Lmax_ref - 5) / 4, tolerance = 1e-10)))
# CR must be CI / RI[n]; canonical Saaty RI table
RI_saaty = c(`1`=0,`2`=0,`3`=0.58,`4`=0.90,`5`=1.12,`6`=1.24,`7`=1.32,
             `8`=1.41,`9`=1.45,`10`=1.49,`11`=1.51,`12`=1.54,`13`=1.56,
             `14`=1.57,`15`=1.59)
check("AHP CR = CI / RI[n]",
      isTRUE(all.equal(unname(r$CR), unname(r$CI / RI_saaty["5"]), tolerance = 1e-10)),
      sprintf("pkg CR=%.6f ref=%.6f", r$CR, r$CI / RI_saaty["5"]))
# perfect consistency: A = w/w^T
w0 = c(0.4, 0.3, 0.2, 0.1)
Apc = outer(w0, w0, "/")
rpc = AHP(Apc)
check("AHP perfect consistency recovers w, CR=0",
      max(abs(rpc$w - w0)) < 1e-8 && abs(rpc$CR) < 1e-8,
      sprintf("maxw=%.2e CR=%.2e", max(abs(rpc$w - w0)), rpc$CR))
# RI[12] exposure: build a 12x12 consistent-ish matrix -> CR = CI/RI[12]
set.seed(7)
A12 = matrix(1, 12, 12); A12[upper.tri(A12)] = runif(66, 1/9, 9)
A12[lower.tri(A12)] = 1 / A12[upper.tri(A12)]; diag(A12) = 1
r12 = AHP(A12)
ev12 = eigen(A12); L12 = Re(ev12$values[which.max(Mod(ev12$values))])
CI12 = (L12 - 12) / 11
cat(sprintf("AHP n=12: CI=%.6f  pkg CR=%.6f -> implied RI=%.6f (canonical 1.54)\n",
            CI12, r12$CR, CI12 / r12$CR))

## ---------------------------------------------------- combine_weights
ws = c(0.4, 0.3, 0.2, 0.1); wo = c(0.25, 0.2, 0.3, 0.25)
check("combine linear",
      isTRUE(all.equal(combine_weights(ws, wo, "linear", 0.6),
                       0.6*ws + 0.4*wo, tolerance = 1e-12)))
check("combine multiplicative (normalized)",
      isTRUE(all.equal(combine_weights(ws, wo, "multiplicative"),
                       (ws*wo)/sum(ws*wo), tolerance = 1e-12)))
check("combine game (geometric mean normalized)",
      isTRUE(all.equal(combine_weights(ws, wo, "game"),
                       sqrt(ws*wo)/sum(sqrt(ws*wo)), tolerance = 1e-12)))
# game_linear: solve crossprod(W) alpha = diag(crossprod(W))
W = cbind(ws, wo); S = crossprod(W)
al = solve(S, diag(S)); al = al/sum(al)
check("combine game_linear",
      isTRUE(all.equal(combine_weights(ws, wo, "game_linear"),
                       as.vector(W %*% al), tolerance = 1e-12)),
      sprintf("pkg=(%s)", paste(round(combine_weights(ws,wo,"game_linear"),5), collapse=",")))

## ------------------------------------------------------- weight_cv
X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
cv_ref = apply(X, 2, \(x) sd(x)/mean(x)); cv_ref = cv_ref/sum(cv_ref)
check("weight_cv", isTRUE(all.equal(unname(weight_cv(X)), as.numeric(cv_ref), tolerance = 1e-12)),
      sprintf("pkg=(%s)", paste(round(weight_cv(X),5), collapse=",")))

## -------------------------------------------------- weight_entropy
Xe = data.frame(x1 = c(3,5,2,7), x2 = c(10,20,15,25))
re = weight_entropy(Xe, c("+","-"))
# independent: rescale to [eps,1-eps], p=x/sum(x), e=-sum(p ln p)/ln(n), w=(1-e)/sum
eps = 0.002
rescale_ref = function(x, type="+", a=eps, b=1-eps){
  rng = range(x); z = (x-rng[1])/(rng[2]-rng[1])
  if(type=="-") z = 1-z
  a + z*(b-a)
}
Xr = cbind(rescale_ref(Xe$x1), rescale_ref(Xe$x2, "-"))
colnames(Xr) = c("x1","x2")
Pm = sweep(Xr, 2, colSums(Xr), "/")
e = -colSums(Pm*log(Pm))/log(nrow(Pm)); w_ref = (1-e)/sum(1-e)
check("weight_entropy w",
      isTRUE(all.equal(unname(re$w), as.numeric(w_ref), tolerance = 1e-10)),
      sprintf("pkg=(%s) ref=(%s)", paste(round(re$w,6),collapse=","), paste(round(w_ref,6),collapse=",")))

## -------------------------------------------------- weight_critic
Xc = data.frame(x1=c(3,5,2,7), x2=c(10,20,15,25))
rc = weight_critic(Xc, c("+","-"))
# independent std-method CRITIC
Xcr = cbind(rescale_ref(Xc$x1, "+", 0, 1), rescale_ref(Xc$x2, "-", 0, 1))
R = rowSums(1 - cor(Xcr)); S = apply(Xcr, 2, sd); C = S*R; w_ref = C/sum(C)
check("weight_critic w (std)",
      isTRUE(all.equal(unname(rc$w), as.numeric(w_ref), tolerance = 1e-10)),
      sprintf("pkg=(%s)", paste(round(rc$w,6),collapse=",")))

## ------------------------------------------------------- linear_sum
Ls = linear_sum(matrix(1:6, 3, 2), c(0.5, 0.5))
check("linear_sum", isTRUE(all.equal(Ls, c((1+4)/2,(2+5)/2,(3+6)/2), tolerance = 1e-12)))

## ------------------------------------------------------- weight_pca
set.seed(1); Xp = data.frame(a=rnorm(30), b=rnorm(30), c=rnorm(30), d=rnorm(30))
rp = weight_pca(Xp, index = rep("+", 4), nfs = 2)
check("weight_pca weights sum to 1", isTRUE(all.equal(sum(rp$w), 1, tolerance = 1e-10)))
check("weight_pca varP = lambda/sum(eigen)", {
  pc = prcomp(scale(Xp), center=FALSE, scale.=FALSE)
  isTRUE(all.equal(unname(rp$lambda), unname((pc$sdev^2)[1:2]), tolerance = 1e-10))
})
## ------------------------------------------------------- weight_fa
set.seed(2); Xf = data.frame(a=rnorm(40), b=rnorm(40), c=rnorm(40))
rf = weight_fa(Xf, nfs = 1)
check("weight_fa weights sum to 1", isTRUE(all.equal(sum(rf$w), 1, tolerance = 1e-10)))
check("weight_fa share = colSums(A^2)/p", {
  A = rf$loading; sh = colSums(A^2)/ncol(Xf)
  isTRUE(all.equal(rf$share, as.vector(sh), tolerance = 1e-10))
})

cat(sprintf("\n=== eval_weights audit: PASS %d / FAIL %d ===\n", P, F))
