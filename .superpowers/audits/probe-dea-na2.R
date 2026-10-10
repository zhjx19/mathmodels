#!/usr/bin/env Rscript
# Verify the user's deterministic super_DEA infeasibility case.
logf = "E:/Github/mathmodels/.superpowers/audits/probe-dea-na2.log"
sink(logf, split = TRUE)
suppressMessages(pkgload::load_all("E:/Github/mathmodels", quiet = TRUE))

dea = data.frame(dmu = LETTERS[1:6],
                 x1 = c(2, 3, 4, 5, 3, 4),
                 x2 = c(3, 4, 5, 6, 4, 5),
                 y1 = c(8, 9, 10, 12, 11, 13),
                 y2 = c(5, 6, 7, 8, 7, 9))
print(dea)

run = function(fn, orient, rts) {
  nw = character(0)
  r = withCallingHandlers(
    tryCatch(fn(dea, inputs = 2:3, outputs = 4:5, orientation = orient, rts = rts),
             error = function(e) paste("ERR:", conditionMessage(e))),
    warning = function(w) { nw <<- c(nw, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.character(r)) { cat("  ERR:", r, "\n"); return(invisible()) }
  nas = names(r$efficiencies)[is.na(r$efficiencies)]
  cat(sprintf("  %-10s %-4s | NA: %-12s | warn: %s\n", fn_name, rts,
              if (length(nas)) paste(nas, collapse = ",") else "(none)",
              if (length(nw)) paste(nw, collapse = " | ") else "(none)"))
}

for (orient in c("io", "oo")) for (rts in c("vrs", "crs")) {
  fn_name = paste0("super_DEA/", orient)
  run(super_DEA, orient, rts)
}

cat("\n--- basic_DEA on the same data (should never fail) ---\n")
for (rts in c("vrs", "crs")) {
  fn_name = "basic_DEA"
  run(basic_DEA, "io", rts)
}

cat("\n--- the user's exact claim ---\n")
r_io = suppressWarnings(super_DEA(dea, 2:3, 4:5, orientation = "io", rts = "vrs"))
r_oo = suppressWarnings(super_DEA(dea, 2:3, 4:5, orientation = "oo", rts = "vrs"))
cat("  io/vrs NA DMU(s):", paste(names(r_io$efficiencies)[is.na(r_io$efficiencies)], collapse = ","), "(claim: F)\n")
cat("  oo/vrs NA DMU(s):", paste(names(r_oo$efficiencies)[is.na(r_oo$efficiencies)], collapse = ","), "(claim: A)\n")
cat("  crs any NA:", any(is.na(suppressWarnings(super_DEA(dea, 2:3, 4:5, orientation="io", rts="crs"))$efficiencies),
                       any(is.na(suppressWarnings(super_DEA(dea, 2:3, 4:5, orientation="oo", rts="crs"))$efficiencies)), "(claim: FALSE)\n")
sink()
cat("\n", readLines(logf, warn = FALSE), sep = "\n")