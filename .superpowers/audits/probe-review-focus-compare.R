# Compare two probe-review-focus-diff.R snapshots (base vs head).
# Usage: Rscript probe-review-focus-compare.R <base_rds> <head_rds>
args = commandArgs(trailingOnly = TRUE)
b = readRDS(args[1]); h = readRDS(args[2])
stopifnot(identical(names(b), names(h)))

same = 0; diff = 0
for (nm in names(b)) {
  x = b[[nm]]; y = h[[nm]]
  ok = identical(x$message, y$message) &&
       identical(x$ok, y$ok) &&
       identical(x$class, y$class) &&
       (if (x$ok) isTRUE(all.equal(x$value, y$value, tolerance = 1e-12)) else TRUE)
  if (ok) { same = same + 1 } else {
    diff = diff + 1
    cat("DIFF:", nm, "\n")
    cat("   base: ok=", x$ok, " class=", paste(x$class, collapse = ","), " msg=", x$message, "\n", sep = "")
    cat("   head: ok=", y$ok, " class=", paste(y$class, collapse = ","), " msg=", y$message, "\n", sep = "")
  }
}
cat(sprintf("\nreview-focus diff: %d cases, %d identical, %d different\n", length(b), same, diff))
