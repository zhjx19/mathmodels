devtools::load_all("E:/Github/mathmodels")
options(warn = -1)
P = 0; F = 0
check = function(n, ok, d="") { if (isTRUE(ok)) {P<<-P+1; cat("PASS:",n,"-",d,"\n")} else {F<<-F+1; cat("FAIL:",n,"-",d,"\n")} }

# SIR closed: S+I+R conserved (no births/deaths)
init = c(S=999, I=1, R=0)
r = ode_sir(init, c(beta=0.3, gamma=0.1), times=seq(0,100,by=1))
total = r$S + r$I + r$R
check("SIR closed conserves N", max(abs(total-1000)) < 1e-3, sprintf("max dev=%.2e", max(abs(total-1000))))
check("SIR I rises then falls", max(r$I) > 1 && r$I[nrow(r)] < max(r$I))

# SEIR closed conserves
init2 = c(S=999, E=0, I=1, R=0)
r2 = ode_seir(init2, c(beta=0.3, gamma=0.1, alpha=0.2), times=seq(0,100,by=1))
t2 = r2$S + r2$E + r2$I + r2$R
check("SEIR closed conserves N", max(abs(t2-1000)) < 1e-3)

# logistic ODE has known solution
rl = ode_logistic(c(N=10), c(r=0.5, K=100), times=seq(0,20,by=1))
K=100; r0=0.5; N0=10
ref = K/(1+(K/N0-1)*exp(-r0*rl$time))
check("logistic matches closed form", isTRUE(all.equal(rl[[2]], ref, tolerance=1e-3)),
      sprintf("maxdiff=%.3e", max(abs(rl[[2]]-ref))))
cat(sprintf("\n=== diffEq audit: PASS %d / FAIL %d ===\n", P, F))
