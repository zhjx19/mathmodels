files = c(list.files("man", pattern = "\\.Rd$", full.names = TRUE),
          list.files("R", pattern = "\\.R$", full.names = TRUE),
          list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE),
          "README.md")
for (f in files) {
  if (!file.exists(f)) next
  lines = readLines(f, encoding = "UTF-8", warn = FALSE)
  for (i in seq_along(lines)) {
    chars = strsplit(lines[i], "")[[1]]
    if (length(chars) == 0) next
    cps = suppressWarnings(vapply(chars, utf8ToInt, integer(1)))
    cps = cps[!is.na(cps) & cps > 127]
    if (length(cps) > 0) {
      cat(sprintf("%s:%d  [%s]  %s\n", f, i,
                  paste(sprintf("U+%04X", unique(cps)), collapse = ","),
                  trimws(lines[i])))
    }
  }
}
cat("scan done\n")
