#!/usr/bin/env Rscript
# Diagnose the R_group1_total values in theil_g2_nest.
logf = "E:/Github/mathmodels/.superpowers/audits/probe-theil.log"
sink(logf, split = TRUE)
suppressMessages(pkgload::load_all("E:/Github/mathmodels", quiet = TRUE))

d = data.frame(
  province = c("A","A","A","A","B","B"),
  city = c("A1","A1","A2","A2","B1","B1"),
  industry = c("Manu","Serv","Manu","Serv","Manu","Serv"),
  y = c(50000, 45000, 60000, 55000, 70000, 65000),
  pop = c(10000, 8000, 15000, 12000, 10000, 8000)
)
r = theil_g2_nest(d, "province", "city", "y", "pop")

cat("ratio:\n"); print(round(r$ratio, 6))
cat("\nthe names of the last block:\n"); print(tail(names(r$ratio), 2))
cat("\nthe values of the last block:\n"); print(tail(r$ratio, 2))

# the .within block is supposed to be: sum over group2 within each group1
cat("\n--- expected per-group1 within contributions ---\n")
# recompute the same way the function does, but correctly grouped
wi = r$ratio[grepl("^[AB]_[A-Z]", names(r$ratio))]
cat("per-(group1,group2) within ratios (Rwi/TT):\n"); print(round(wi, 6))
cat("  sum =", round(sum(wi), 6), " (= ratio['Tw'] )\n")
cat("  group A alone =", round(sum(wi[grepl('^A_', names(wi))]), 6), "\n")
cat("  group B alone =", round(sum(wi[grepl('^B_', names(wi))]), 6), "\n")
cat("\n>> .within block returned:\n"); print(round(tail(r$ratio, 2), 6))
cat(">> each returned value equals the GLOBAL within total, not the group's own:\n")
cat("   group A .within == total within ?",
    isTRUE(all.equal(unname(tail(r$ratio, 2)[1]), sum(wi))), "\n")
sink()
cat(readLines(logf, warn = FALSE), sep = "\n")