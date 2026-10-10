#!/usr/bin/env Rscript
# Run one testthat filter against an arbitrary package dir (for a pristine baseline).
a = commandArgs(trailingOnly = TRUE)
pkgdir = a[[1]]; filter = a[[2]]; logf = a[[3]]
sink(logf)
options(warn = 1)
cat("pkg:", pkgdir, "| filter:", filter, "\n")
r = devtools::test(pkgdir, filter = filter, stop_on_failure = FALSE,
                   reporter = "summary")
sink()
d = as.data.frame(r)
tot = vapply(c("failed", "error", "warning", "skipped", "passed"),
             function(k) if (k %in% names(d)) sum(as.numeric(d[[k]])) else 0,
             numeric(1))
cat(sprintf("\nSUMMARY failed=%d error=%d warning=%d skipped=%d passed=%d\n",
            tot[["failed"]], tot[["error"]], tot[["warning"]],
            tot[["skipped"]], tot[["passed"]]),
    file = logf, append = TRUE)