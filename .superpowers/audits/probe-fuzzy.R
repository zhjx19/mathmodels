# Independent audit probe: eval_fuzzy.R
devtools::load_all("E:/Github/mathmodels")
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}

## -------- membership functions: textbook shape checks --------
# triangular
x = seq(0, 10, 0.5)
y = mf_tri(x, c(3, 6, 8))
check("mf_tri peak=1 at b, 0 at a,c",
      y[x==6] == 1 && y[x==3] == 0 && y[x==8] == 0)
check("mf_tri linear rises", isTRUE(all.equal(
  y[x>3 & x<6], (x[x>3 & x<6]-3)/3, tolerance=1e-12)))
check("mf_tri linear falls", isTRUE(all.equal(
  y[x>6 & x<8], (8-x[x>6 & x<8])/2, tolerance=1e-12)))
check("mf_tri range [0,1]", all(y >= 0 & y <= 1))

# trapezoidal
yt = mf_trap(x, c(1, 5, 7, 8))
check("mf_trap plateau=1", all(yt[x>5 & x<7] == 1))
check("mf_trap shoulders", yt[x==1]==0 && yt[x==8]==0 && abs(yt[x==5]-1)<1e-12)
# one-sided shoulder
ysh = mf_trap(x, c(-Inf, -Inf, 3, 5))
check("mf_trap left shoulder", all(ysh[x<3]==1) && ysh[x==5]==0)

# gaussian
yg = mf_gauss(x, c(2, 5))
check("mf_gauss peak=1 at c", yg[x==5] == 1)
check("mf_gauss formula", isTRUE(all.equal(
  yg, exp(-(x-5)^2/(2*2^2)), tolerance=1e-12)))

# generalized bell
ygb = mf_gbell(x, c(2, 4, 6))
check("mf_gbell peak=1 at c", abs(ygb[x==6]-1) < 1e-12)
check("mf_gbell formula", isTRUE(all.equal(
  ygb, 1/(1+abs((x-6)/2)^(2*4)), tolerance=1e-12)))

# two-param gaussian
yg2 = mf_gauss2(seq(0,10,0.5), c(1, 3, 3, 4))
check("mf_gauss2 plateau between c1,c2", all(yg2[seq(0,10,0.5)>=3 & seq(0,10,0.5)<=4] == 1))
check("mf_gauss2 left decay", abs(yg2[seq(0,10,0.5)==0] - exp(-9/2)) < 1e-12)

# sigmoid family
check("sigmoid(0, a=1,b=0)=0.5", abs(sigmoid(0,1,0)-0.5) < 1e-12)
sd1 = mf_sigmoid(x, c(2,4)); check("mf_sigmoid at b=0.5", abs(mf_sigmoid(4,c(2,4))-0.5)<1e-12)
sd2 = mf_dsigmoid(x, c(5,2,5,7))
check("mf_dsigmoid = sig - sig", isTRUE(all.equal(
  sd2, sigmoid(x,5,2)-sigmoid(x,5,7), tolerance=1e-12)))
sp = mf_psigmoid(x, c(2,3,-5,8))
check("mf_psigmoid = sig * sig", isTRUE(all.equal(
  sp, sigmoid(x,2,3)*sigmoid(x,-5,8), tolerance=1e-12)))

# z / pi / s shapes
yz = mf_z(x, c(3,7))
check("mf_z: 1 left, 0 right, 0.5 midpoint",
      yz[x<=3][1]==1 && yz[x>7][1]==0 && abs(yz[x==5]-0.5)<1e-12)
yp = mf_pi(x, c(1,4,5,10))
check("mf_pi plateau and tangents",
      all(yp[x>=4 & x<=5]==1) && yp[x==1]==0 && abs(yp[x==10])<1e-12)
ys = mf_s(x, c(1,8))
check("mf_s: 0 left, 1 right, 0.5 midpoint",
      ys[x<=1][1]==0 && ys[x>=8][1]==1 && abs(ys[x==4.5]-0.5)<1e-12)
check("all mf in [0,1]", all(yz>=0&yz<=1) && all(yp>=0&yp<=1) && all(ys>=0&ys<=1))

## -------- build_mf_tri / gauss / sigmoid + compute_mf --------
th = c(0.05, 0.15, 0.25, 0.5)
mfs_tri = compute_mf_funs(th, "tri")
check("compute_mf_funs tri returns 4 funcs", length(mfs_tri) == 4)
mm = compute_mf(0.07, th, "tri")
# level1: trap(-Inf,-Inf,0.05,0.15) -> decay 1->0
check("compute_mf tri level1", abs(mm[1] - (0.15-0.07)/(0.15-0.05)) < 1e-12)
# level2: tri(0.05,0.15,0.25) -> rise
check("compute_mf tri level2", abs(mm[2] - (0.07-0.05)/(0.15-0.05)) < 1e-12)
check("compute_mf len = n knots", length(mm) == 4)
mfs_g = compute_mf_funs(th, "gauss", sigma=0.05)
mg = compute_mf(0.07, th, "gauss", sigma=0.05)
check("compute_mf gauss level2 = gauss centered 0.15",
      abs(mg[2] - exp(-(0.07-0.15)^2/(2*0.05^2))) < 1e-12)
mfs_s = compute_mf_funs(th, "sigmoid", slope=20)
ms = compute_mf(0.07, th, "sigmoid", slope=20)
check("compute_mf sigmoid first level decreasing at th1",
      ms[1] < 0.5, sprintf("lvl1=%.4f", ms[1]))
check("compute_mf all in [0,1]", all(mm>=0&mm<=1) && all(mg>=0&mg<=1) && all(ms>=0&ms<=1))

## -------- fuzzy_eval operators vs hand formulas --------
w = c(0.3, 0.3, 0.3, 0.1)
R = matrix(c(0.8,0.7,0.6,0.7, 0.1,0.2,0.2,0.1, 0.1,0.1,0.2,0.2), nrow=3, byrow=TRUE)
b3 = fuzzy_eval(w, R, 3)
B3 = rowSums(sweep(R, 2, w, "*"))
check("fuzzy_eval type 3 = weighted sum normalized",
      isTRUE(all.equal(b3, B3/sum(B3), tolerance=1e-12)),
      sprintf("pkg=(%s)", paste(round(b3,5),collapse=",")))
b1 = fuzzy_eval(w, R, 1)
B1 = apply(R, 1, \(r) max(pmin(w, r)))
check("fuzzy_eval type 1 = max(min(w,r))",
      isTRUE(all.equal(b1, B1/sum(B1), tolerance=1e-12)))
b2 = fuzzy_eval(w, R, 2)
B2 = apply(R, 1, \(r) max(w*r))
check("fuzzy_eval type 2 = max(w*r)", isTRUE(all.equal(b2, B2/sum(B2), tolerance=1e-12)))
b4 = fuzzy_eval(w, R, 4)
B4 = apply(R, 1, \(r) min(1, sum(pmin(w, r))))
check("fuzzy_eval type 4 = min(1,sum(min))", isTRUE(all.equal(b4, B4/sum(B4), tolerance=1e-12)))
b5 = fuzzy_eval(w, R, 5)
B5 = apply(R, 1, \(r) sum(pmin(w, r/sum(r))))
check("fuzzy_eval type 5 = sum(min(w,r/sum r))", isTRUE(all.equal(b5, B5/sum(B5), tolerance=1e-12)))
check("fuzzy_eval B sums to 1", isTRUE(all.equal(sum(b3), 1, tolerance=1e-12)))

## -------- defuzzify --------
mu = c(0.318, 0.351, 0.203, 0.128); sc = c(30,60,75,90)
check("defuzzify weighted_average", abs(defuzzify(mu, sc, "weighted_average") - sum(mu*sc)) < 1e-12)
check("defuzzify centroid", abs(defuzzify(mu, sc, "centroid") - sum(mu*sc)/sum(mu)) < 1e-12)
check("defuzzify max_membership", defuzzify(mu, sc, "max_membership") == 60)
# centroid with unnormalized mu
mu2 = c(2, 1); check("centroid normalizes unnormalized mu",
  abs(defuzzify(mu2, c(10,20), "centroid") - (20+20)/3) < 1e-12)

cat(sprintf("\n=== eval_fuzzy audit: PASS %d / FAIL %d ===\n", P, F))
