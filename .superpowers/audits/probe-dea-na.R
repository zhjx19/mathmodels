#!/usr/bin/env Rscript
# Try to force solve_radial_lp to return status != 0 (=> silent NA).
options(warn = 1)
logf = "E:/Github/mathmodels/.superpowers/audits/probe-dea-na.log"
sink(logf, split = TRUE)
suppressMessages(pkgload::load_all("E:/Github/mathmodels", quiet = TRUE))

chk = function(label, df, inputs, outputs, fn = super_DEA, ...) {
  nw = 0
  r = withCallingHandlers(
    tryCatch(fn(df, inputs = inputs, outputs = outputs, ...),
             error = function(e) paste("ERR:", conditionMessage(e))),
    warning = function(wl) { nw <<- nw + 1
                             cat("   [warn]", conditionMessage(wl), "\n")
                             invokeRestart("muffleWarning") })
  if (is.character(r)) { cat(sprintf("  %-38s %s\n", label, r)); return(invisible()) }
  nas = names(r$theta)[is.na(r$theta)]
  cat(sprintf("  %-38s NA: %s | warnings: %d\n", label,
              if (length(nas)) paste(nas, collapse = ",") else "(none)", nw))
}

cat("\n== super_DEA: look for an infeasible LP ==\n")
# DMU1 strictly dominates (min input, max output) -> its exclusion may be infeasible
d1 = data.frame(DMU = paste0("D", 1:5), x = c(1, 5, 5, 5, 5), y = c(10, 5, 5, 5, 5))
chk("dominating D1, vrs", d1, 2, 3, super_DEA, rts = "vrs")
chk("dominating D1, crs", d1, 2, 3, super_DEA, rts = "crs")

# two identical extreme DMUs
d2 = data.frame(DMU = paste0("D", 1:4), x = c(1, 1, 9, 9), y = c(10, 10, 1, 1))
chk("two extreme DMUs, vrs", d2, 2, 3, super_DEA, rts = "vrs")

# zero input (LP bounds force lambda >= 0; x0 = 0)
d3 = data.frame(DMU = paste0("D", 1:4), x = c(0, 3, 6, 9), y = c(9, 6, 3, 0))
chk("zero input DMU, vrs", d3, 2, 3, super_DEA, rts = "vrs")
chk("zero input DMU, crs", d3, 2, 3, super_DEA, rts = "crs")

d4 = data.frame(DMU = paste0("D", 1:4), x = c(0, 3, 6, 9), y = c(0, 3, 6, 9))
chk("zero input + zero output, vrs", d4, 2, 3, super_DEA, rts = "vrs")

cat("\n== basic_DEA: same silent-NA pattern? ==\n")
chk("basic_DEA dom D1 vrs", d1, 2, 3, basic_DEA, rts = "vrs")
chk("basic_DEA zero-input vrs", d3, 2, 3, basic_DEA, rts = "vrs")

cat("\n== many random draws, both models ==\n")
set.seed(42)
hit_s = hit_b = 0
for (k in 1:2000) {
  n = sample(3:8, 1)
  df = data.frame(DMU = paste0("D", 1:n),
                  x1 = runif(n, 0, 10), x2 = runif(n, 0, 10),
                  y1 = runif(n, 0, 10))
  for (rts in c("vrs", "crs")) {
    s = suppressWarnings(tryCatch(super_DEA(df, 2:3, 4, rts = rts),
                                  error = function(e) NULL))
    if (!is.null(s) && anyNA(s$theta)) hit_s = hit_s + 1
    b = suppressWarnings(tryCatch(basic_DEA(df, 2:3, 4, rts = rts),
                                  error = function(e) NULL))
    if (!is.null(b) && anyNA(b$theta)) hit_b = hit_b + 1
  }
}
cat("  super_DEA draws with NA:", hit_s, " basic_DEA draws with NA:", hit_b, "\n")

cat("\n=== done ===\n")
sink()
cat(readLines(logf, warn = FALSE), sep = "\n")