#!/usr/bin/env Rscript
# Reproduce every claim in mathmodels-0.0.14-修复补丁建议.md against source.
options(warn = 1)
logf = "E:/Github/mathmodels/.superpowers/audits/verify-patch-claims.log"
sink(logf, split = TRUE)

pkg = "E:/Github/mathmodels"
suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
cat("mathmodels loaded from:", pkg, "\n\n")

hr   = function(x) cat("\n", strrep("-", 6), x, strrep("-", 6), "\n")
safe = function(label, expr) {
  withCallingHandlers(
    tryCatch(expr, error = function(e)
      cat("  !! ERROR:", conditionMessage(e), "\n")),
    warning = function(wl) { cat("  [warning]", conditionMessage(wl), "\n")
                             invokeRestart("muffleWarning") })
}

# --- 1. fuzzy_eval(type = 5) with a zero-membership grade -------------------
hr("1. fuzzy_eval type 5, grade 1 all-zero membership")
safe("1", {
  w  = c(0.2, 0.1, 0.15, 0.3, 0.25)
  R5 = matrix(c(0, 0, 0, 0, 0,
                0.3, 0.4, 0.5, 0.2, 0.6,
                0.7, 0.6, 0.5, 0.8, 0.4), nrow = 3, byrow = TRUE)
  cat("type = 3 ->", paste(round(fuzzy_eval(w, R5, 3), 4), collapse = ", "), "\n")
  t5 = tryCatch(fuzzy_eval(w, R5, 5),
                error = function(e) paste("ERR:", conditionMessage(e)))
  cat("type = 5 ->", if (is.character(t5)) t5 else paste(round(t5, 4), collapse = ", "), "\n")
  r0 = sum(R5[1, ])
  Bj = suppressWarnings(sum(apply(rbind(w, R5[1, ] / r0), 2, min)))
  cat("  r0 =", r0, "| B[1] =", Bj, "| is.nan(B[1]) =", is.nan(Bj), "\n")
  cat("  `NaN == 0` ->", suppressWarnings(NaN == 0), "\n")
  cat("  >> package's own 'Normalization failed' message reached:",
      is.character(t5) && grepl("Normalization failed", t5), "\n")
  t0 = tryCatch(fuzzy_eval(c(0.5, 0.5), matrix(0, 2, 2), 5),
                error = function(e) paste("ERR:", conditionMessage(e)))
  cat("  all-zero R ->", t0, "\n")
})

# --- 2. theil_g2_nest ratio structure ---------------------------------------
hr("2. theil_g2_nest(): ratio length / duplicate names")
safe("2", {
  xl = "D:/WorkBuddy/mathmodels技能/mathmodels-book/data/income.xlsx"
  if (!file.exists(xl)) { cat("  SKIP: income.xlsx not found\n") } else {
    d = readxl::read_excel(xl)
    ycol = if ("平均工资" %in% names(d)) "平均工资" else "y"
    pcol = if ("就业人数" %in% names(d)) "就业人数" else "pop"
    g1 = "区域"; g2 = "地区"
    cat("  group1 =", g1, "| group2 =", g2, "| y =", ycol, "| pop =", pcol, "\n")
    r = theil_g2_nest(d, g1, g2, ycol, pcol)
    cat("  length(theil) =", length(r$theil),
        "| length(ratio) =", length(r$ratio),
        "| diff =", length(r$ratio) - length(r$theil), "\n")
    dup = unique(names(r$ratio)[duplicated(names(r$ratio))])
    cat("  duplicated names in ratio:", paste(dup, collapse = ", "), "\n")
    cat("  ratio[['东北']] =", r$ratio[["东北"]], "\n")
    cat("  >> lengths equal:", length(r$theil) == length(r$ratio), "\n")
    tt = tryCatch({ tibble::tibble(a = unname(r$theil), b = unname(r$ratio)); "ok" },
                  error = function(e) conditionMessage(e))
    cat("  >> tibble(theil, ratio) side by side:", tt, "\n")
  }
})

# --- 3. mf_dsigmoid range ---------------------------------------------------
hr("3. mf_dsigmoid() range vs documented [0,1]")
safe("3", {
  x = seq(0, 10, by = 0.01)
  for (p in list(c(5, 2, 5, 7), c(5, 7, 5, 2), c(1, 1.5, 9, 1.5))) {
    v = mf_dsigmoid(x, p)
    cat(sprintf("  params %-20s range = [%8.4f, %8.4f]  any negative: %s\n",
                paste(p, collapse = ","), min(v), max(v), any(v < 0)))
  }
})

# --- 4. markov_chain dimnames ----------------------------------------------
hr("4. markov_chain(): unnamed dimnames")
safe("4", {
  S  = c("A", "B", "A", "C", "A", "B", "C", "B", "A", "C")
  mc = markov_chain(S, "A", 3)
  P  = mc$trans_mat
  cat("  class(P):", paste(class(P), collapse = "/"), "\n")
  cat("  names(dimnames(P)):", paste0("'", names(dimnames(P)), "'", collapse = ","),
      "| all empty:", all(!nzchar(names(dimnames(P)))), "\n")
  tb = tryCatch(tibble::as_tibble(P, rownames = "x"),
                error = function(e) paste("ERR:", conditionMessage(e)))
  cat("  as_tibble(P, rownames='x') ->",
      if (is.character(tb)) tb else "ok (no error)", "\n")
})

# --- 5. weight_entropy upper-bound guard -----------------------------------
hr("5. weight_entropy(): no upper-bound warning on un-normalised NA column")
safe("5", {
  X  = data.frame(a = 1:10, b = 1:10)   # b marked NA but NOT normalised to [0,1]
  nw = 0
  wr = withCallingHandlers(weight_entropy(X, index = c("+", NA)),
        warning = function(wl) { nw <<- nw + 1
                                 cat("  WARNING:", conditionMessage(wl), "\n")
                                 invokeRestart("muffleWarning") })
  cat("  number of warnings:", nw, "| range(s) = [",
      round(min(wr$s), 2), ",", round(max(wr$s), 2), "]\n")
  cat("  >> silent for values > 1:", nw == 0, "\n")
})

# --- 6. super_DEA silent NA -------------------------------------------------
hr("6. super_DEA(): silent NA when LP fails")
safe("6", {
  hit = NULL
  for (k in 1:400) {
    n  = 5
    dd = data.frame(DMU = paste0("D", 1:n),
                    x = runif(n, 0.1, 10), y = runif(n, 0.1, 10))
    nw = 0
    sr = withCallingHandlers(
      super_DEA(dd, inputs = 2, outputs = 3, rts = "crs"),
      warning = function(wl) { nw <<- nw + 1; invokeRestart("muffleWarning") })
    if (anyNA(sr$theta)) { hit = list(k = k, dmu = names(sr$theta)[is.na(sr$theta)],
                                      warn = nw); break }
  }
  if (is.null(hit))
    cat("  no NA in 400 random trials; code path still lacks an else-branch\n")
  else
    cat("  NA at trial", hit$k, "for DMU(s):", paste(hit$dmu, collapse = ", "),
        "| warnings raised:", hit$warn, "\n")
})

cat("\n=== done ===\n")
sink()
cat(readLines(logf, warn = FALSE), sep = "\n")