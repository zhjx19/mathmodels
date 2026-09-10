# Independent audit probe: eval_algorithms + eval_system + eval_regional
devtools::load_all("E:/Github/mathmodels")
options(warn = -1)
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}

## ================= TOPSIS =================
X = data.frame(X1 = c(2,5,3,4), X2 = c(8,1,6,2), X3 = c(5,5,3,7))
w = c(0.5, 0.3, 0.2); idx = c("+","-","+")
# independent TOPSIS
B = sweep(as.matrix(X), 2, sqrt(colSums(X^2)), "/")
C = sweep(B, 2, w, "*")
mx = apply(C, 2, max); mn = apply(C, 2, min)
Cstar = ifelse(idx=="+", mx, mn); C0 = ifelse(idx=="+", mn, mx)
Ss = sqrt(rowSums((C - matrix(Cstar, nrow(C), ncol(C), byrow=TRUE))^2))
S0 = sqrt(rowSums((C - matrix(C0, nrow(C), ncol(C), byrow=TRUE))^2))
ref = S0/(S0+Ss); names(ref) = paste0("Sample", 1:4)
check("topsis vs independent implementation",
      isTRUE(all.equal(unname(topsis(X, w, idx)), unname(ref), tolerance = 1e-12)),
      sprintf("pkg=(%s) ref=(%s)", paste(round(topsis(X,w,idx),5),collapse=","), paste(round(ref,5),collapse=",")))
check("topsis all-positive default", isTRUE(all.equal(
  unname(topsis(X, w)), unname(topsis(X, w, rep("+",3))), tolerance=1e-12)))

## ================= grey_corr =================
refg = c(0.9, 0.8, 0.7)
cmp = data.frame(x1=c(0.9,0.7,0.8), x2=c(0.8,0.9,0.7), x3=c(0.7,0.8,0.9))
g = grey_corr(refg, cmp, rho=0.5)
# independent GRA
p = nrow(cmp); wv = rep(1/p, p)
T = sweep(as.matrix(cmp), 1, refg, "-")
min2 = min(abs(T)); max2 = max(abs(T))
eta = (min2 + 0.5*max2)/(abs(T) + 0.5*max2)
ref_g = as.vector(wv %*% eta)
check("grey_corr vs independent GRA",
      isTRUE(all.equal(unname(g), ref_g, tolerance = 1e-12)),
      sprintf("pkg=(%s)", paste(round(g,5),collapse=",")))

## ================= rank_sum_ratio =================
rsd = data.frame(ID=c("A","B","C"), X1=c(10,20,15), X2=c(5,10,8))
rr = rank_sum_ratio(rsd, c(0.4,0.6), "int")
# independent weighted RSR (int ranks, ascending -> larger better)
rk = sapply(rsd[-1], rank)          # 3 x 2
rsr_ref = as.vector(rk %*% c(0.4,0.6)) / (1*3)
check("rank_sum_ratio RSR values",
      isTRUE(all.equal(sort(rr$resultTable$RSR), sort(rsr_ref), tolerance=1e-12)),
      sprintf("pkg RSR=(%s) ref=(%s)", paste(round(sort(rr$resultTable$RSR),5),collapse=","),
              paste(round(sort(rsr_ref),5),collapse=",")))
check("rank_sum_ratio returns lm and tables",
      inherits(rr$reg, "lm") && all(c("resultTable","rankTable") %in% names(rr)))

## ================= coupling_degree (standard) =================
df = data.frame(s1=c(0.2,0.8,0.5), s2=c(0.4,0.6,0.5), s3=c(0.3,0.9,0.5))
cd = coupling_degree(df, id_cols=NULL)
# independent standard coupling C = p*(prod x)^(1/p)/sum(x); T = weighted mean; D=sqrt(C*T)
p3 = 3; w3 = rep(1/3,3)
Cref = apply(df, 1, \(x) p3*prod(x)^(1/p3)/sum(x))
Tref = as.vector(as.matrix(df) %*% w3)
Dref = sqrt(Cref*Tref)
check("coupling_degree standard C", isTRUE(all.equal(cd$coupling, Cref, tolerance=1e-12)))
check("coupling_degree standard T,D", isTRUE(all.equal(cd$coord, Tref, tolerance=1e-12)) &&
      isTRUE(all.equal(cd$coupling_coord, Dref, tolerance=1e-12)))
check("coupling_degree: equal subsystems -> C=1", {
  de = data.frame(a=c(0.5,0.5), b=c(0.5,0.5), c=c(0.5,0.5))
  all(abs(coupling_degree(de)$coupling - 1) < 1e-12)
})
check("coupling_degree adjusted equal -> C=1", {
  de = data.frame(a=0.5, b=0.5, c=0.5)
  abs(coupling_degree(de, type="adjusted")$coupling - 1) < 1e-12
})
check("coupling_degree adjusted in [0,1]", {
  ca = coupling_degree(df, type="adjusted")$coupling
  all(ca >= -1e-12 & ca <= 1+1e-12)
})

## ================= obstacle_degree =================
od = obstacle_degree(df, c(0.4,0.3,0.3), id_cols=NULL)
Oref = sweep(1-as.matrix(df), 2, c(0.4,0.3,0.3), "*")
check("obstacle_degree O=(1-X)*w", isTRUE(all.equal(as.matrix(od), Oref, tolerance=1e-12)))

## ================= LQ =================
lqd = data.frame(region=c("National","A","B"), total=c(10000,4000,6000),
                 ind1=c(2000,1000,1000), ind2=c(6000,2000,4000))
lq = LQ(lqd, region, total, dplyr::starts_with("ind"))
# independent: (region share of ind)/(national share of ind)
lq_ref = (lqd$ind1[-1]/lqd$total[-1]) / (lqd$ind1[1]/lqd$total[1])
check("LQ ind1", isTRUE(all.equal(lq$ind1, lq_ref, tolerance=1e-12)),
      sprintf("pkg=(%s) ref=(%s)", paste(round(lq$ind1,5),collapse=","), paste(round(lq_ref,5),collapse=",")))

## ================= HHI =================
check("HHI raw", isTRUE(all.equal(HHI(c(50,30,20)), (0.5^2+0.3^2+0.2^2), tolerance=1e-12)))
check("HHI scaled", isTRUE(all.equal(HHI(c(50,30,20), TRUE),
        ((0.5^2+0.3^2+0.2^2) - 1/3)/(1-1/3), tolerance=1e-12)))
check("HHI single =1", HHI(5) == 1)

## ================= EG =================
eg = data.frame(
  region=c("R1","R1","R1","R1","R2","R2","R3","R3","R1","R2","R2","R3"),
  industry=c("A","A","A","A","A","A","A","A","B","B","B","B"),
  emp=c(250,200,150,100,20,15,10,5,50,200,150,50))
reg = EG(eg, region, industry, emp)
# independent EG for each industry
eg_ref = sapply(sort(unique(eg$industry)), function(ind) {
  di = eg[eg$industry==ind,]
  xr = tapply(eg$emp, eg$region, sum); xr = xr/sum(xr)   # region share of total emp
  sr = tapply(di$emp, di$region, sum); sr = sr/sum(sr)   # region share of industry emp
  # align regions
  regs = sort(unique(eg$region))
  x = as.numeric(xr[regs]); s = as.numeric(sr[regs]); s[is.na(s)] = 0
  h = sum((di$emp/sum(di$emp))^2)
  if (h==1) 0 else (sum((s-x)^2)/(1-sum(x^2)) - h)/(1-h)
})
check("EG vs independent formula",
      isTRUE(all.equal(unname(reg$EG), unname(eg_ref), tolerance = 1e-10)),
      sprintf("pkg=(%s) ref=(%s)", paste(round(reg$EG,6),collapse=","), paste(round(eg_ref,6),collapse=",")))

cat(sprintf("\n=== algorithms/system/regional audit: PASS %d / FAIL %d ===\n", P, F))
