# Differential characterisation probe for the Phase A equivalence refactor.
#
#   capture     : run the grid on the CURRENT build and snapshot to
#                 probe-validation-baseline.rds
#   compare     : re-run the grid and assert everything matches the snapshot,
#                 except that the error *class* may gain "rlang_error".
#   attribution : assert the reported call site of representative errors names
#                 the EXPORTED function. Only meaningful AFTER migration: sites
#                 that used `stop(..., call. = FALSE)` had a NULL call before.
#
# Run via: Rscript probe-validation-equivalence.R capture|compare|attribution
devtools::load_all("E:/Github/mathmodels", quiet = TRUE)

BASE = "E:/Github/mathmodels/.superpowers/sdd/eval-validation-plan/probe-validation-baseline.rds"

## --- input grid -------------------------------------------------------------
df4  = data.frame(a = c(8, 7, 6, 9), b = c(150, 180, 200, 160), c = c(60, 80, 100, 70))
df4n = df4; df4n$a[2] = NA
dff  = data.frame(a = letters[1:3], b = 1:3)          # non-numeric column
M4   = as.matrix(df4)
one  = data.frame(a = 1, b = 2, c = 3)                # 1 row
empty = data.frame(a = numeric(0), b = numeric(0), c = numeric(0))
w3   = c(0.3, 0.4, 0.3)
i3   = c("+", "+", "+")
i3na = c("+", NA, "+")
i3bad = c("+", "x", "+")
w_inf = c(0.4, 0.3, Inf)                              # Review Focus #1
w_nan = c(0.4, 0.3, NaN)                              # Review Focus #1

cases = list(
  # matrix/data.frame guard: df-or-matrix sites
  "topsis:str"        = quote(topsis("a")),
  "topsis:df"         = quote(topsis(df4, w3, i3)),
  "topsis:M"          = quote(topsis(M4, w3, i3)),
  "topsis:nonnum"     = quote(topsis(dff, w3, i3)),            # Review Focus #2
  "topsis:1row"       = quote(topsis(one, w3, i3)),            # Review Focus #4
  "topsis:empty"      = quote(topsis(empty, w3, i3)),          # Review Focus #5
  "topsis:w_len"      = quote(topsis(df4, c(0.5), i3)),
  "topsis:w_na"       = quote(topsis(df4, c(0.5, NA), i3)),
  "topsis:w_inf"      = quote(topsis(df4, w_inf, i3)),
  "topsis:w_nan"      = quote(topsis(df4, w_nan, i3)),
  "topsis:i_len"      = quote(topsis(df4, w3, c("+", "+"))),
  "topsis:i_bad"      = quote(topsis(df4, w3, i3bad)),
  "topsis:i_na"       = quote(topsis(df4, w3, i3na)),
  # grey_corr_topsis
  "gct:df"            = quote(grey_corr_topsis(df4, w3, i3)),
  "gct:M"             = quote(grey_corr_topsis(M4, w3, i3)),
  "gct:i_len"         = quote(grey_corr_topsis(df4, w3, c("+", "+"))),
  "gct:i_nulllen"     = quote(grey_corr_topsis(df4, w3, NULL)),
  # grey_corr / rank_sum_ratio
  "gc:str"            = quote(grey_corr(c(1, 2, 3), "a")),
  "gc:len"            = quote(grey_corr(c(1, 2), M4)),
  "rsr:str"           = quote(rank_sum_ratio(1:5)),
  # weights family
  "critic:str"        = quote(weight_critic("a")),
  "critic:df"         = quote(weight_critic(df4, i3na)),
  "critic:M"          = quote(weight_critic(M4, i3na)),
  "critic:nonnum"     = quote(weight_critic(dff, i3)),         # Review Focus #2
  "critic:1row"       = quote(weight_critic(one, i3)),         # Review Focus #4
  "critic:i_len"      = quote(weight_critic(df4, c("+", "+"))),
  "critic:i_bad"      = quote(weight_critic(df4, i3bad)),
  "cv:str"            = quote(weight_cv("a")),
  "cv:df"             = quote(weight_cv(df4)),
  "cv:M"              = quote(weight_cv(as.matrix(abs(df4)))),
  "entropy:str"       = quote(weight_entropy("a")),
  "entropy:df"        = quote(weight_entropy(df4, i3na)),
  "entropy:M"         = quote(weight_entropy(abs(M4), i3)),
  "fa:str"            = quote(weight_fa("a")),
  "fa:df"             = quote(weight_fa(df4, i3, nfs = 1)),
  "fa:M"              = quote(weight_fa(M4, i3, nfs = 1)),
  "pca:str"           = quote(weight_pca("a")),
  "pca:df"            = quote(weight_pca(df4, i3na, nfs = 1)),
  "pca:M"             = quote(weight_pca(M4, i3, nfs = 1)),
  "pca:i_len"         = quote(weight_pca(df4, c("+", "+"), nfs = 1)),
  "linear_sum:str"    = quote(linear_sum("a", c(1, 2))),
  "linear_sum:len"    = quote(linear_sum(df4, c(0.4, 0.3, 0.3))),
  # system
  "coupling:str"      = quote(coupling_degree("a")),
  "coupling:df"       = quote(coupling_degree(df4, w3)),
  "coupling:w_len"    = quote(coupling_degree(df4, c(0.5, 0.5))),
  "coupling:w_sum"    = quote(coupling_degree(df4, c(0.5, 0.4, 0.4))),
  "obstacle:str"      = quote(obstacle_degree("a")),
  "obstacle:w_len"    = quote(obstacle_degree(df4, c(0.5, 0.5))),
  # fuzzy
  "fuzzy:w_num"       = quote(fuzzy_eval("a", matrix(1, 3, 3), 1)),
  "fuzzy:w_len"       = quote(fuzzy_eval(c(0.3, 0.3), matrix(1, 3, 3), 1)),
  "fuzzy:R_df"        = quote(fuzzy_eval(w3, as.data.frame(matrix(0.5, 3, 3)), 1)),
  # dea / inequality / regional
  "dea:str"           = quote(basic_DEA("a", inputs = 1, outputs = 2)),
  "malmquist:df"      = quote(malmquist(subset(mtcars, select = c(mpg, cyl, hp)), period = "cyl", inputs = 1, outputs = 2)),
  "theil_g:str"       = quote(theil_g(1:5, "g", "y", "p")),
  "theil0_g:str"      = quote(theil0_g(1:5, "g", "y")),
  "LQ:str"            = quote(LQ("a", region, total, .cols = x)),
  "EG:str"            = quote(EG("a", region, industry, y))
)

run_case = function(expr) {
  err = NULL
  val = withCallingHandlers(
    tryCatch(eval(expr, envir = globalenv()), error = function(e) { err <<- e; NULL }),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (!is.null(err)) {
    list(ok = FALSE, value = NULL,
         class = class(err),
         message = conditionMessage(err),
         call = if (is.null(conditionCall(err))) NA_character_ else paste(deparse(conditionCall(err)), collapse = " "))
  } else {
    list(ok = TRUE,
         value = if (is.numeric(val)) signif(as.numeric(val), 12) else val,
         class = class(val), message = NA_character_, call = NA_character_)
  }
}

mode = commandArgs(trailingOnly = TRUE)[1]
results = lapply(cases, run_case)
names(results) = names(cases)

if (identical(mode, "capture")) {
  saveRDS(results, BASE)
  cat("captured", length(results), "cases ->", BASE, "\n")
} else if (identical(mode, "compare")) {
  base = readRDS(BASE)
  bad = 0
  for (nm in names(base)) {
    b = base[[nm]]; n = results[[nm]]
    if (!identical(b$message, n$message)) {
      cat("MISMATCH message:", nm, "\n  base:", b$message, "\n  new :", n$message, "\n"); bad = bad + 1
    } else if (b$ok != n$ok) {
      cat("MISMATCH ok:", nm, "\n"); bad = bad + 1
    } else if (b$ok && !isTRUE(all.equal(b$value, n$value, tolerance = 1e-12))) {
      cat("MISMATCH value:", nm, "\n"); bad = bad + 1
    } else if (!b$ok && !all(b$class %in% c("error", "simpleError")) &&
               !all(n$class %in% c(b$class, "rlang_error"))) {
      cat("MISMATCH class:", nm, "\n  base:", paste(b$class, collapse = ","),
          "\n  new :", paste(n$class, collapse = ","), "\n"); bad = bad + 1
    }
  }
  cat(sprintf("compare: %d cases, %d mismatches\n", length(base), bad))
  if (bad > 0) quit(status = 1)
} else if (identical(mode, "attribution")) {
  # Post-migration only: the reported call must name the EXPORTED function.
  expect_fn = c(
    "topsis:str" = "topsis", "topsis:i_bad" = "topsis", "topsis:w_na" = "topsis",
    "gct:i_len" = "grey_corr_topsis", "critic:str" = "weight_critic",
    "critic:i_len" = "weight_critic", "pca:i_len" = "weight_pca",
    "entropy:str" = "weight_entropy", "cv:str" = "weight_cv",
    "coupling:str" = "coupling_degree", "obstacle:str" = "obstacle_degree",
    "linear_sum:str" = "linear_sum", "theil_g:str" = "theil_g",
    "theil0_g:str" = "theil0_g", "LQ:str" = "LQ", "EG:str" = "EG",
    # `extract_dea_data` is DEA's internal choke point and already owned this
    # error before migration (baseline records the same call) -- unchanged.
    "dea:str" = "extract_dea_data",
    "fa:str" = "weight_fa", "rsr:str" = "rank_sum_ratio"
  )
  bad = 0
  for (nm in names(expect_fn)) {
    cl = results[[nm]]$call
    if (is.na(cl) || !grepl(expect_fn[[nm]], cl, fixed = TRUE)) {
      cat("ATTRIBUTION MISS:", nm, "expected", expect_fn[[nm]], "got:", cl, "\n"); bad = bad + 1
    }
  }
  cat(sprintf("attribution: %d cases, %d misses\n", length(expect_fn), bad))
  if (bad > 0) quit(status = 1)
} else {
  stop("usage: probe-validation-equivalence.R capture|compare|attribution")
}
