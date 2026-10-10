#!/usr/bin/env Rscript
# Full R CMD check for mathmodels, output captured to a log.
pkg  = "E:/Github/mathmodels"
logf = "E:/Github/mathmodels/.superpowers/audits/r-cmd-check.log"

sink(logf)
cat("=== devtools::check() ===\n")
cat("started:", format(Sys.time()), "\n")
cat("R:", R.version.string, "| devtools:", as.character(packageVersion("devtools")), "\n\n")

res = tryCatch(
  devtools::check(pkg, document = FALSE, cran = FALSE, error_on = "never"),
  error = function(e) { cat("\n*** check() itself errored:", conditionMessage(e), "\n"); NULL }
)

cat("\n\n=== returned object ===\n")
if (!is.null(res)) {
  print(res)
  cat("\nerrors:", length(res$errors), "| warnings:", length(res$warnings),
      "| notes:", length(res$notes), "\n")
  if (length(res$errors))   { cat("\n--- ERRORS ---\n");   print(res$errors) }
  if (length(res$warnings)) { cat("\n--- WARNINGS ---\n"); print(res$warnings) }
  if (length(res$notes))    { cat("\n--- NOTES ---\n");    print(res$notes) }
} else {
  cat("no result object\n")
}
cat("\nfinished:", format(Sys.time()), "\n")
sink()