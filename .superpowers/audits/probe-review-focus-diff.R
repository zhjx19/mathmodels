# Base-vs-head differential probe for the plan's Review Focus input classes.
#
# The Phase A plan claimed its 57-case grid covered all five Review Focus
# classes. It did not: the two `*:nonnum` cases were mis-sized (a 2-column data
# frame with a length-3 `index`), so they hit the index-length guard and never
# reached the non-numeric column. This probe compares a BASE revision against
# HEAD for correctly-sized inputs, so the claim is either substantiated or
# falsified on evidence.
#
# Usage:
#   Rscript probe-review-focus-diff.R <pkg_root> <out_rds>
#   Rscript probe-review-focus-compare.R <base_rds> <head_rds>
devtools_load = function(root) devtools::load_all(root, quiet = TRUE)

args = commandArgs(trailingOnly = TRUE)
root = args[1]; out = args[2]
devtools_load(root)

df4  = data.frame(a = c(8, 7, 6, 9), b = c(150, 180, 200, 160), c = c(60, 80, 100, 70))
M4   = as.matrix(df4)
# I1: a VALID matrix that happens to contain a zero in one column
Mzero = matrix(c(0, 7, 6, 9, 150, 180, 200, 160, 60, 80, 100, 70), 4, 3)
# Review Focus #2: non-numeric column, CORRECTLY sized (2 cols -> index length 2)
dff2 = data.frame(a = letters[1:3], b = 1:3)
i2   = c("+", "+")
one  = data.frame(a = 1, b = 2, c = 3)
zerorow = data.frame(a = numeric(0), b = numeric(0), c = numeric(0))

cases = list(
  # --- RF#1: Inf / NaN in w (not NA) ---
  "rf1:topsis_w_inf"    = quote(topsis(df4, c(0.4, 0.3, Inf), c("+", "+", "+"))),
  "rf1:topsis_w_nan"    = quote(topsis(df4, c(0.4, 0.3, NaN), c("+", "+", "+"))),
  "rf1:gct_w_inf"       = quote(grey_corr_topsis(df4, c(0.4, 0.3, Inf), c("+", "+", "+"))),
  "rf1:critic_bad_idx"  = quote(weight_critic(df4, c(0.4, 0.3, Inf))),
  # --- RF#2: non-numeric column, correctly sized so the guard is passed ---
  "rf2:topsis_nonnum"   = quote(topsis(dff2, c(0.5, 0.5), i2)),
  "rf2:critic_nonnum"   = quote(weight_critic(dff2, i2)),
  "rf2:entropy_nonnum"  = quote(weight_entropy(dff2, i2)),
  "rf2:pca_nonnum"      = quote(weight_pca(dff2, i2)),
  "rf2:fa_nonnum"       = quote(weight_fa(dff2, i2, nfs = 1)),
  "rf2:linear_nonnum"   = quote(linear_sum(dff2, c(0.5, 0.5))),
  "rf2:gct_nonnum"      = quote(grey_corr_topsis(dff2, c(0.5, 0.5), i2)),
  # --- RF#3: index mixing NA with + / - ---
  "rf3:critic_i_na"     = quote(weight_critic(df4, c("+", NA, "+"))),
  "rf3:entropy_i_na"    = quote(weight_entropy(df4, c("+", NA, "+"))),
  "rf3:pca_i_na"        = quote(weight_pca(df4, c("+", NA, "+"))),
  "rf3:gct_i_na"        = quote(grey_corr_topsis(df4, c(0.4, 0.3, 0.3), c("+", NA, "+"))),
  # --- RF#4: 1-row / 1-column ---
  "rf4:topsis_1row"     = quote(topsis(one, c(0.4, 0.3, 0.3), c("+", "+", "+"))),
  "rf4:critic_1row"     = quote(weight_critic(one, c("+", "+", "+"))),
  "rf4:pca_1row"        = quote(weight_pca(one, c("+", "+", "+"))),
  "rf4:entropy_1col"    = quote(weight_entropy(data.frame(a = c(1, 2, 3)))),
  # --- RF#5: 0-row / 0-col ---
  "rf5:topsis_0row"     = quote(topsis(zerorow, c(0.4, 0.3, 0.3), c("+", "+", "+"))),
  "rf5:critic_0row"     = quote(weight_critic(zerorow, c("+", "+", "+"))),
  "rf5:linear_0row"     = quote(linear_sum(zerorow, c(0.4, 0.3, 0.3))),
  # --- I1: matrix with a zero vs the identical data frame ---
  "i1:topsis_matrix"    = quote(topsis(Mzero, c(0.3, 0.4, 0.3), c("+", "+", "+"))),
  "i1:topsis_dataframe" = quote(topsis(as.data.frame(Mzero), c(0.3, 0.4, 0.3), c("+", "+", "+")))
)

run_case = function(expr) {
  err = NULL
  val = withCallingHandlers(
    tryCatch(eval(expr, envir = globalenv()), error = function(e) { err <<- e; NULL }),
    warning = function(w) invokeRestart("muffleWarning")
  )
  if (!is.null(err)) {
    list(ok = FALSE, value = NULL, class = class(err),
         message = conditionMessage(err),
         call = if (is.null(conditionCall(err))) NA_character_ else paste(deparse(conditionCall(err)), collapse = " "))
  } else {
    list(ok = TRUE,
         value = if (is.numeric(val)) signif(as.numeric(val), 12) else val,
         class = class(val), message = NA_character_, call = NA_character_)
  }
}

saveRDS(lapply(cases, run_case), out)
cat("wrote", length(cases), "cases from", root, "->", out, "\n")
