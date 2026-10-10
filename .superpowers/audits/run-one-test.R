#!/usr/bin/env Rscript
# Run one testthat filter. Usage: Rscript run-one-test.R <filter> <logfile>
a = commandArgs(trailingOnly = TRUE)
filter = a[[1]]
logf = a[[2]]

sink(logf)
options(warn = 1)
r = devtools::test(filter = filter, stop_on_failure = FALSE, reporter = "summary")
sink()

d = as.data.frame(r)
tot = vapply(c("failed", "error", "warning", "skipped", "passed"),
             function(k) if (k %in% names(d)) sum(as.numeric(d[[k]])) else 0,
             numeric(1))
cat(sprintf("\nSUMMARY failed=%d error=%d warning=%d skipped=%d passed=%d\n",
            tot[["failed"]], tot[["error"]], tot[["warning"]],
            tot[["skipped"]], tot[["passed"]]),
    file = logf, append = TRUE)