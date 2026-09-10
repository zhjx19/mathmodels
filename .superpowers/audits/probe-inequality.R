# Independent audit probe: eval_inequality.R
devtools::load_all("E:/Github/mathmodels")
P = 0; F = 0
check = function(name, ok, detail = "") {
  if (isTRUE(ok)) { P <<- P + 1; cat("PASS:", name, "-", detail, "\n") }
  else { F <<- F + 1; cat("FAIL:", name, "-", detail, "\n") }
}

## trapz vs pracma
x = c(0, 1, 2, 3.5); y = c(0, 2, 1, 4)
check("trapz vs pracma::trapz",
      isTRUE(all.equal(trapz(x, y), pracma::trapz(x, y), tolerance = 1e-12)),
      sprintf("pkg=%.10f pracma=%.10f", trapz(x, y), pracma::trapz(x, y)))

## gini0 vs references
set.seed(1); inc = c(10, 20, 30, 40, 100)
check("gini0 vs ineq::Gini",
      isTRUE(all.equal(gini0(inc), ineq::Gini(inc), tolerance = 1e-10)),
      sprintf("pkg=%.8f ineq=%.8f", gini0(inc), ineq::Gini(inc)))
check("gini0 vs DescTools::Gini (sample-corrected convention)",
      isTRUE(all.equal(gini0(inc), DescTools::Gini(inc) * (length(inc)-1)/length(inc),
                       tolerance = 1e-10)),
      sprintf("pkg=%.8f Desc*n/(n-1) ratio=%.8f", gini0(inc),
              DescTools::Gini(inc) * (length(inc)-1)/length(inc)))
# (DescTools::Gini defaults to the n/(n-1) sample-corrected estimator;
#  ineq::Gini and gini0 use the population formula - convention difference)
# random vector
inc2 = c(1.3, 5.7, 2.2, 9.9, 4.1, 7.3, 0.5)
check("gini0 vs ineq::Gini (random)",
      isTRUE(all.equal(gini0(inc2), ineq::Gini(inc2), tolerance = 1e-10)))

## theil0 vs references (Theil T)
check("theil0 vs ineq::Theil (T)",
      isTRUE(all.equal(theil0(inc), ineq::Theil(inc), tolerance = 1e-10)),
      sprintf("pkg=%.8f ineq=%.8f", theil0(inc), ineq::Theil(inc)))
check("theil0 vs ineq::Theil (random)",
      isTRUE(all.equal(theil0(inc2), ineq::Theil(inc2), tolerance = 1e-10)))

## singleton groups: grouped == individual
check("gini singleton == gini0",
      isTRUE(all.equal(gini(inc, rep(1, length(inc))), gini0(inc), tolerance = 1e-12)),
      sprintf("gini=%.8f gini0=%.8f", gini(inc, rep(1, 5)), gini0(inc)))
check("theil singleton == theil0",
      isTRUE(all.equal(theil(inc, rep(1, length(inc))), theil0(inc), tolerance = 1e-12)),
      sprintf("theil=%.8f theil0=%.8f", theil(inc, rep(1, 5)), theil0(inc)))

## gini grouped: independent Lorenz-curve implementation
gini_ref = function(y, pop) {
  o = order(y / pop); y = y[o]; pop = pop[o]
  x = c(0, cumsum(pop)/sum(pop)); z = c(0, cumsum(y)/sum(pop*y/sum(pop)))
  z = c(0, cumsum(y)/sum(y))
  1 - sum((x[-1]-x[-length(x)])*(z[-1]+z[-length(z)]))
}
incg = c(10, 20, 30, 40, 100); popg = c(100, 150, 200, 250, 300)
check("gini grouped vs independent Lorenz",
      isTRUE(all.equal(gini(incg, popg), gini_ref(incg, popg), tolerance = 1e-12)),
      sprintf("pkg=%.8f ref=%.8f", gini(incg, popg), gini_ref(incg, popg)))

## theil grouped: independent implementation
theil_ref = function(y, pop) {
  Yb = sum(y*pop)/sum(pop); Yr = y/Yb
  sum(pop/sum(pop) * Yr * ifelse(Yr > 0, log(Yr), 0))
}
check("theil grouped vs independent",
      isTRUE(all.equal(theil(incg, popg), theil_ref(incg, popg), tolerance = 1e-12)))

## ------- decomposition identity tests -------
## theil0_g: TT must equal theil0 on the full micro data
data = data.frame(g = c("A","A",rep("B",10),rep("A",6)),
                  y = c(10,10,rep(8,4),rep(6,6),rep(4,4),2,2))
r0 = theil0_g(data, "g", "y")
check("theil0_g: T == theil0(micro)",
      isTRUE(all.equal(r0$theil["T"], c(T = theil0(data$y))[1], tolerance = 1e-10)),
      sprintf("pkg T=%.8f direct=%.8f", r0$theil["T"], theil0(data$y)))
check("theil0_g: T == Tb + Tw",
      isTRUE(all.equal(unname(r0$theil["T"]),
                       unname(r0$theil["Tb"] + r0$theil["Tw"]), tolerance = 1e-12)),
      sprintf("T=%.8f Tb+Tw=%.8f", r0$theil["T"], r0$theil["Tb"] + r0$theil["Tw"]))

## theil_g: TT must equal theil on the grouped-average data
data2 = data |> dplyr::count(g, y, name = "pop")
rg = theil_g(data2, "g", "y", "pop")
check("theil_g: T == theil(grouped)",
      isTRUE(all.equal(unname(rg$theil["T"]), theil(data2$y, data2$pop), tolerance = 1e-10)),
      sprintf("pkg T=%.8f direct=%.8f", rg$theil["T"], theil(data2$y, data2$pop)))
check("theil_g: T == Tb + Tw",
      isTRUE(all.equal(unname(rg$theil["T"]),
                       unname(rg$theil["Tb"] + rg$theil["Tw"]), tolerance = 1e-12)))

## theil_g2_cross: TT must equal theil on raw rows; T == Tb + Tw_b + Tw_w
data3 = data.frame(
  industry = c("A","A","A","A","B","B","B","B"),
  area = c("East","East","West","West","East","East","West","West"),
  province = c("Shanghai","Beijing","Sichuan","Yunnan","Shanghai","Beijing","Sichuan","Yunnan"),
  avg_wage = c(500,400,80,60,300,250,50,40),
  emp_num = c(100,80,90,70,120,100,80,60))
rc = theil_g2_cross(data3, "industry", "area", "avg_wage", "emp_num")
check("theil_g2_cross: TT == theil(raw rows)",
      isTRUE(all.equal(unname(rc$theil["T"]), theil(data3$avg_wage, data3$emp_num), tolerance = 1e-10)),
      sprintf("pkg TT=%.8f direct=%.8f", rc$theil["T"], theil(data3$avg_wage, data3$emp_num)))
check("theil_g2_cross: TT == Tb + Tw_b + Tw_w",
      isTRUE(all.equal(unname(rc$theil["T"]),
                       unname(rc$theil["Tb"] + rc$theil["Tw_b"] + rc$theil["Tw_w"]), tolerance = 1e-10)),
      sprintf("TT=%.8f sum=%.8f", rc$theil["T"], rc$theil["Tb"] + rc$theil["Tw_b"] + rc$theil["Tw_w"]))

## theil_g2_nest: TT must equal theil on raw rows; T == Tb_g1 + Tb_g2 + Tw
data4 = data.frame(
  province = c("A","A","A","A","B","B"),
  city = c("A1","A1","A2","A2","B1","B1"),
  industry = c("Manu","Serv","Manu","Serv","Manu","Serv"),
  y = c(50000,45000,60000,55000,70000,65000),
  pop = c(10000,8000,15000,12000,10000,8000))
rn = theil_g2_nest(data4, "province", "city", "y", "pop")
check("theil_g2_nest: TT == theil(raw rows)",
      isTRUE(all.equal(unname(rn$theil["T"]), theil(data4$y, data4$pop), tolerance = 1e-10)),
      sprintf("pkg TT=%.8f direct=%.8f", rn$theil["T"], theil(data4$y, data4$pop)))
check("theil_g2_nest: TT == Tb_g1 + Tb_g2 + Tw",
      isTRUE(all.equal(unname(rn$theil["T"]),
                       unname(rn$theil["Tb_g1"] + rn$theil["Tb_g2"] + rn$theil["Tw"]), tolerance = 1e-10)),
      sprintf("TT=%.8f sum=%.8f", rn$theil["T"],
              rn$theil["Tb_g1"] + rn$theil["Tb_g2"] + rn$theil["Tw"]))

## edges
check("gini0 rejects negative", inherits(tryCatch(gini0(c(1,-2,3)), error = \(e) e), "error"))
check("theil0 handles zero (0*log0 = 0)", is.finite(theil0(c(0,1,2))))
check("gini handles zero-income group", is.finite(gini(c(0,1,2), c(1,1,1))))

cat(sprintf("\n=== eval_inequality audit: PASS %d / FAIL %d ===\n", P, F))
