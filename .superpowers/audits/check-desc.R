#!/usr/bin/env Rscript
# CRAN DESCRIPTION sanity checks per cran-extrachecks skill.
pkg = "E:/Github/mathmodels"
dcf = read.dcf(file.path(pkg, "DESCRIPTION"))
field = function(n) if (n %in% colnames(dcf)) dcf[1, n] else NA_character_

ok = TRUE
say = function(pass, label, detail = "") {
  if (!pass) ok <<- FALSE
  cat(sprintf("[%s] %s%s\n", if (pass) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0(" -> ", detail) else ""))
}

# --- parses at all ---
say(TRUE, "read.dcf() parses DESCRIPTION")

# --- authoritative: the same checker R CMD check uses ---
rr = tools:::.check_package_description(file.path(pkg, "DESCRIPTION"))
say(length(rr) == 0, "tools:::.check_package_description: 0 issues",
    sprintf("%d issue(s)", length(rr)))

title = field("Title")
desc = field("Description")
authors = field("Authors@R")
lic = field("License")

# --- Title ---
say(!grepl("(^|\\s)[Aa] (Toolkit|Tools) for", title), "Title: no 'A Toolkit/Tools for'")
say(!grepl("(?i)\\bfor R\\b", title), "Title: no trailing 'for R'")
say(nchar(title) <= 65, "Title: <= 65 chars", sprintf("%d chars", nchar(title)))
tc = tools::toTitleCase(gsub("[A-Za-z][A-Za-z']*", "\\0", title))
say(TRUE, "Title: suggested title case", tc)

# --- Description ---
# strip continuation lines -> single logical string
lines = strsplit(desc, "\n")[[1]]
flat = trimws(paste(lines, collapse = " "))
flat = gsub("[[:space:]]+", " ", flat)

first = sub("[[:space:]].*$", "", flat)
say(!grepl("^(This package|This|Package) ", first) &&
      !grepl("^[A-Za-z.]+ (provides|offers|implements) ", first) &&
      !grepl("^Functions? for ", first),
    "Description: no forbidden opening", first)

sents = gregexpr("[.](\\s|$)", flat)[[1]]
nsent = length(sents)
say(nsent >= 3 && nsent <= 5, "Description: 3-5 sentences",
    sprintf("%d sentences", nsent))

say(!grepl(",\\s*,|,\\s*\\.", flat), "Description: no comma splices")
say(!grepl(",\\s*,", desc), "Description: no ',,' artifact")
say(!grepl("[.?!]\\s*,\\s*(and|or)\\s*[a-z]", desc), "Description: no '. , and' artifact")
say(!grepl('"', flat), "Description: no double quotes")

# acronyms: capitalised tokens of length >= 2 that are not real words
acronyms = regmatches(flat, gregexpr("\\b[A-Z]{2,}\\b", flat))[[1]]
known = c("DE", "Word", "LaTeX")
unexp = setdiff(unique(acronyms), known)
say(length(unexp) == 0, "Description: no unexplained acronyms",
    if (length(unexp)) paste(unexp, collapse = ", ") else "none")

# software names must be single-quoted; function names must NOT be
say(!grepl("'[A-Za-z_.]+\\(\\)'", desc), "Description: function names not quoted")

# any bare standalone token 'R' must appear as 'R'
bare_r = regmatches(flat, gregexpr("(?<!')\\bR\\b", flat, perl = TRUE))[[1]]
say(length(bare_r) == 0, "Description: no unquoted bare 'R'",
    if (length(bare_r)) paste(bare_r, collapse = ", ") else "none")

# --- Authors / license ---
say(grepl("cph", authors), "Authors@R: has cph role")
say(grepl("cre", authors), "Authors@R: has cre role")
say(nzchar(lic), "License: present", lic)

# --- URLs ---
urls = unlist(strsplit(paste(field("URL"), field("BugReports"), field("Repository")),
                         "\\s*,\\s*"))
urls = urls[nzchar(urls)]
say(all(grepl("^https://", urls)), "All URLs use https",
    paste(urls, collapse = " | "))

cat("\n--- flattened Description ---\n")
cat(flat, "\n\n")
cat(if (ok) "RESULT: all checks passed\n" else "RESULT: FAILURES above\n")
quit(status = if (ok) 0 else 1)