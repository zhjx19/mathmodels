#!/usr/bin/env Rscript
logf = "E:/Github/mathmodels/.superpowers/audits/full-suite.log"
sink(logf)
options(warn = 1)
r = devtools::test("E:/Github/mathmodels", stop_on_failure = FALSE, reporter = "summary")
sink()
d = as.data.frame(r)
tot = vapply(c("failed", "error", "warning", "skipped", "passed", "nb"),
             function(k) if (k %in% names(d)) sum(as.numeric(d[[k]])) else 0,
             numeric(1))
cat(sprintf("files=%d  failed=%d error=%d warning=%d skipped=%d passed=%d\n",
            nrow(d), tot[["failed"]], tot[["error"]], tot[["warning"]],
            tot[["skipped"]], tot[["passed"]]),
    file = logf, append = TRUE)