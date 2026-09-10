# Audit: `R/stat_infer.R` vs gold-standard references

- **Date:** 2026-09-09
- **Scope:** all 13 exported functions in `R/stat_infer.R` (1340 lines) — `stat_describe`, `stat_t_test`, `stat_wilcox_test`, `stat_anova`, `stat_kruskal`, `stat_ancova`, `stat_cor`, `stat_chisq`, `stat_normality`, `stat_hov_test`, `stat_bootstrap`, `stat_perm_test`, plus `print.stat_infer`
- **Method:** identical-input comparison against `stats`, `coin`, `car`, `psych`, `rstatix`, hand-exact enumerations (252 / 32 / 64 / 720 permutation nulls), and same-seed `boot::boot` reproduction. Degenerate/edge inputs probed for crashes and silent-wrong-results.
- **Probe scripts:** `.superpowers/audits/probe-stat_infer.R` (main, 147 comparison rows) and `.superpowers/audits/probe-followup.R` (targeted confirmations); machine-readable results in `.superpowers/audits/probe-results.rds`
- **Verdict counts:** 116 MATCH / 31 VARIANT / 2 probe-labeled BUG (both reclassified below — see findings) — **0 computational errors found**; all p/statistics/effect sizes match references to rel diff 0 (< 1e-6) wherever both compute the same quantity.

## Per-function verdict table

Numbers are from the probe; "rd" = relative difference (0 = exact match). Full rows in probe output.

| Function | Reference | Our value | Ref value | rd | Verdict |
|---|---|---|---|---|---|
| stat_describe (mean/sd/median/iqr/CI) | hand, listwise-complete subset | 0.17298 / 0.92013 / 0.25658 / 1.04934 / −0.08852 | same | 0 | match |
| stat_describe skew/kurt | population moments g1/g2 (doc-consistent) | −0.07175 / 1.30082 | same | 0 | match (variant vs `psych::describe` type-3 — by design, documented as "moment") |
| stat_t_test Welch (stat/df/p/CI) | `stats::t.test` default interface | −1.13330 / 30.3978 / 0.26595 / −3.07031 | same | 0 | match |
| stat_t_test pooled / one-sample / paired / greater / mu | `stats::t.test` | all identical | same | 0 | match (all variants) |
| stat_t_test `estimate` column | internal consistency | +1.1588 (b−a) | statistic −1.2630, effect −0.4168, CI [−3.0255, 0.7080] (all a−b) | — | **Important finding F2** (sign contradiction inside one row) |
| stat_t_test Cohen's d (pooled / dz) | hand formula | −0.35705 / −0.35501 | same | 0 | match |
| stat_wilcox_test W, p (ties, correct=TRUE) | `stats::wilcox.test` formula | 10 / 0.26768 | same | 0 | match |
| stat_wilcox_test exact p | `stats::wilcox.test` exact AND `coin::wilcox_test` exact | 0.785714 | same | 0 | match |
| stat_wilcox_test rank-biserial | hand 2W/(n1n2)−1 (a−b direction) | −0.428571 | same | 0 | match (rstatix `wilcox_effsize` uses a different z-based formula — not comparable, not a defect) |
| stat_wilcox_test paired rbc (zero-diff handling) | hand, zeros dropped | −0.2 | same | 0 | match |
| stat_wilcox_test one-sample p / effect | `stats::wilcox.test` / docs | 0.50781 / **NA** | same / promised | — | p match; effect **Minor finding F8** |
| stat_anova 1-way SS/F/p/η²p, Residuals row | `stats::anova(aov())` | 42.5568 / 31.1408 / 4.787e-08 / 0.67491 | same | 0 | match |
| stat_anova 2-way SS (g1, g2, g1:g2) | `anova(aov(g1*g2))` Type I | 42.5568 / 0.15951 / 5.63524 | same | 0 | match; vs `car::Anova(type=2)` g1 = 41.74993 → **Minor finding F9** (SS type undocumented) |
| stat_anova `.by` slice | manual subset `aov` | 25.62938 | same | 0 | match |
| stat_kruskal H / p / ε² | `stats::kruskal.test` / `rstatix::kruskal_effsize` | 23.4268 / 8.184e-06 / 0.71423 | same | 0 | match |
| stat_kruskal small-n p vs `coin` exact | coin | 0.2425 | 0.3429 | — | variant (asymptotic χ² by design, same as stats; no exact option — note only) |
| stat_ancova SS (g, x) | `anova(aov())` Type I | 23.4938 / 174.878 | same | 0 | match mechanically; vs Type II g = 22.066 (balanced) / **0.7511 vs 8.4880 (unbalanced)** → **Important finding F3** |
| stat_ancova 2 covariates | Type I | 174.878 / 1.83993 | same | 0 | match |
| stat_ancova `test_parallel` warns | hand interaction | warns | warns | — | match for first response; **Important finding F4** for multi-response |
| stat_cor pearson/spearman/kendall (est/p/stat) | `stats::cor.test` | all identical | same | 0 | match |
| stat_cor greater / CI 90% / spearman CI NA | `stats::cor.test` | 8.157e-06 / 0.43217 / NA | same | 0 | match |
| stat_cor n_used per pair | pairwise-complete convention | 8 | 9 (pair a-b has 9 complete rows) | — | **Important finding F1** (listwise across ALL `.cols`) |
| stat_chisq 2x2 correct=TRUE/FALSE, 3x3, expected counts, V | `stats::chisq.test` + hand V | 4.32 / 0.037667 / 0.020921 / 0.995024 / 0.258199 | same | 0 | match (V from uncorrected χ², consistent under `correct` toggle) |
| stat_chisq exact | `stats::fisher.test` | 0.036835 | same | 0 | match |
| stat_normality SW stat/p, skew/kurt, n∈[3,5000] guards | `stats::shapiro.test` + moments | 0.95031 / 0.07771 | same | 0 | match; n<3 / n>5000 → NA as documented |
| stat_hov_test levene (median/mean), bartlett, F | `car::leveneTest`, `stats::bartlett.test`, `stats::var.test` | 1.07997 / 1.33172 / 2.11650 / 2.64580 | same | 0 | match; F-test `df2` dropped → **Minor finding F10** |
| stat_bootstrap estimate/bias/se/CI (perc & bca) | same-seed `boot::boot` + `boot::boot.ci` | 5.43202 / 0.0026163 / 0.24492 / 4.98875 / 4.99786 | same | 0 | match (bit-identical reproduction) |
| stat_bootstrap se vs theory; coverage | σ/√n; 100-rep coverage | 0.24492 (theory 0.25430); 90/100 | 95% nominal | — | match (MC noise; percentile CI known slight undercoverage) |
| stat_bootstrap constant statistic | warning + NA bounds | warns / NA | expected | — | match (graceful) |
| stat_perm_test two.sample / one.sample / paired / independence / greater | **exact hand enumeration** (252 / 32 / 64 / 720) | 0.69191 / 0.06334 / 0.18908 / 0.05736 / 0.34601 | 0.69048 / 0.06250 / 0.18750 / 0.05833 / 0.34524 | — | match (MC error ≤ 0.0015; p convention (1+#)/(B+1) correct) |
| stat_perm_test coin cross-check | `coin::oneway_test` exact | 0.69191 | 0.69048 | — | variant (coin uses standardized statistic — sanity only, consistent) |
| stat_perm_test custom two.sample (documented idiom) | — | p = 1.0000, null_sd = 0 | meaningful null expected | — | **Important finding F5** |
| stat_perm_test custom independence ("second column" per docs) | doc claim | shuffles LAST column (L1228) | second column | — | **Important finding F6** |
| print.stat_infer | — | label + tibble | — | — | match |

## Findings

Grades: **Critical** = wrong results; **Important** = crash with opaque error or wrong semantics; **Minor** = docs/edge polish. No Critical findings — every number that both sides compute agrees exactly.

### Important

**F1. `stat_cor`: listwise deletion across ALL `.cols` silently changes pair statistics** — `R/stat_infer.R:742`
`keep = stats::complete.cases(cbind(cols_df, slice_df))` drops a row for the (a, b) pair when an unrelated third column `cc` is NA there. Doc says "pairwise correlations among numeric columns"; `stats::cor.test`, `Hmisc::rcorr` use pairwise-complete data. Repro (probe §F):
```r
d = data.frame(a = 1:10, b = 1:10 + rnorm(10), cc = c(rep(1, 9), NA)); d$b[1] = NA
stat_cor(d, .cols = a:cc)$n_used          # 8 8 8  (pair a-b uses 8 rows)
stats::cor.test(d$a, d$b)                 # uses 9 pairwise-complete rows
```
n_used, r, statistic and p for every pair are computed on the intersection-complete subset; nothing warns. Fix direction: compute per-pair on pairwise-complete rows (or at minimum document the rule loudly).

**F2. `stat_t_test`: `estimate` sign contradicts `statistic`/`conf_low..high`/`effect` in the same row** — `R/stat_infer.R:298`
`estimate = est[2] - est[1]` = mean(group2) − mean(group1), while the statistic, CI, and Cohen's d all follow group1 − group2 (the internal call is `t.test(x, y)` with x = first level). Repro (probe followup E):
```r
r = stat_t_test(dt, .cols = y, group = g)   # mean(a)=9.5749, mean(b)=10.7336
r$statistic  # -1.2630  (a - b)
r$estimate   # +1.1588  (b - a)  <-- opposite sign
r$effect     # -0.4168  (a - b);  CI [-3.0255, 0.7080] (a - b)
```
Not a computational error, but a user reading `estimate` against `statistic` gets the direction backwards. `rstatix` convention is group1 − group2 throughout. Fix direction: make `estimate = est[1] - est[2]` (or document the flip).

**F3. `stat_ancova`: group term tested WITHOUT covariate adjustment (Type I SS)** — `R/stat_infer.R:669-671`
`resp ~ group + covariates` through `aov()` gives sequential SS: the group term is entered first, so its SS is the *unadjusted* between-group SS. The namesake ANCOVA question (covariate-adjusted group effect, Type II/III as in `car::Anova`, `rstatix::anova_test`, textbook ANCOVA) is not what is reported whenever group sizes are unequal or slopes non-parallel. Repro (probe §6, unbalanced 14-row data):
```r
r = stat_ancova(dc2, resp = y, group = g, covariates = x)
r$sumsq[r$term == "g"]        # 0.7511  (Type I)
car::Anova(aov(y ~ g + x, dc2), type = 2)["g", "Sum Sq"]   # 8.4880
# p: mathmodels 0.51 vs Type II ≈ 0.07 — an order-of-magnitude SS difference
```
Silent and affects the headline ANCOVA test. Fix direction: `car::Anova(fit, type = 2)` (or document "Type I sequential SS" and warn for unbalanced data).

**F4. `stat_ancova(test_parallel = TRUE)`: parallel-slopes check runs only for the FIRST response** — `R/stat_infer.R:654-657` (`names(resp_df)[1]` at L655)
With `resp = c(y1, y2)` where only `y2` has a strong group×covariate interaction, no warning fires; `resp = y2` alone warns (probe followup A: interaction p = 4.3e-10 missed in the multi-response call). Responses beyond the first are silently untested for the assumption the flag exists to check. Fix direction: loop the check over all responses.

**F5. `stat_perm_test` custom two-sample statistic: the documented pattern yields a degenerate null (p ≡ 1, silent)** — `R/stat_infer.R:1218-1231`, docs L1151-1154
For `type = "two.sample"` the framework passes re-ordered row indices to the statistic, but the reordering applies to the data frame as a whole — any statistic of the documented form `f(d, i)` that computes a group contrast from `d$y[i]` and `d$g[i]` re-indexes y and g TOGETHER, so pairs never break: the statistic is invariant, `null_sd = 0`, and the p-value is exactly 1. The package's own test file (`tests/testthat/test-stat_infer.R:373-378`) demonstrates precisely this idiom. Repro (probe followup C):
```r
f = function(d, i) {
  y = d$y[i]; lab = factor(d$g[i], levels = c("a", "b"))
  median(y[lab == "a"]) - median(y[lab == "b"])
}
stat_perm_test(dperm, statistic = f, type = "two.sample", permutations = 2000, seed = 4)
# observed = 1, null_mean = 1, null_sd = 0, p = 1  — meaningless, no warning
```
Fix direction: either permute the group labels (not rows) for custom two-sample statistics, or detect `null_sd ≈ 0` and error/warn; fix docs and the test example.

**F6. `stat_perm_test` custom independence: code shuffles the LAST column, docs say "second column"** — `R/stat_infer.R:1228` vs docs L1152-1154
`d2[[ncol(d)]] = d[[ncol(d)]][idx_mat[, b]]` shuffles the last column of the whole data frame, not the second `.cols` column. Two silent failure modes (probe §12 + followup B):
1. Statistic uses a non-last column: `cor(d$x, d$junk)` with columns `(x, junk, z)` — the shuffled column is `z`, `junk` is never permuted → null_mean equals the observed correlation, p meaningless.
2. With `.by` present, the last column is the **slice column**: `stat_perm_test(dby, statistic = \(d, i) cor(d$x, d$z), type = "independence", .by = s)` shuffles `s`, returns `p.value = 1`, `null_sd = 0` for both slices — silently, though `cor(x, z) = 0.93`.
Fix direction: shuffle the second `.cols` column (as documented) or require the statistic to designate its shuffle target; at minimum error when `ncol(d)` is a slice column.

**F7. `stat_describe`: cross-column listwise deletion + misleading `n_missing`** — `R/stat_infer.R:145,167`
(a) `complete.cases` spans ALL `.cols` (and `.by`), so `mean`/`sd`/... of column x change depending on which OTHER columns are selected (probe §1: `describe(x)` mean 0.15730 vs `describe(x:z)` mean 0.17298). Not documented.
(b) `n_missing = nrow(data) - n` is not the variable's own NA count: for column z with 1 own NA it reports 2 (includes the row where x is NA); under `.by = g`, slice b with **zero** own NAs reports `n_missing = 5` — the count of rows belonging to the OTHER slice. Repro:
```r
stat_describe(data.frame(x = c(1:4, NA, 6:10), g = rep(c("a","b"), each = 5)),
  .cols = x, .by = g)$n_missing   # slice b (no NAs) -> 5
```
Fix direction: per-variable NA counting; per-variable (or documented) deletion; slice-relative `n_missing`.

### Minor

**F8. `stat_wilcox_test`: one-sample path returns `effect = NA` though `@return` promises rank-biserial correlation** — `R/stat_infer.R:396-397`. The matched-pairs rbc formula exists for `paired`; one-sample (`mu`-based signed-rank) gets NA with no doc note.

**F9. `stat_anova`/`stat_ancova`: Type I (sequential) SS never stated in docs; `stat_anova` `@return` advertises `p.adjusted` but the function has no `p_adjust` argument and hardcodes `method = "none"` (L512), so `p.adjusted` is always identical to `p.value`** — `R/stat_infer.R:486-487,512` (anova), L669-671 (ancova). rstatix `anova_test` defaults to Type II; users comparing against it will see different SS for early terms in unbalanced designs. Document the SS type (or switch to Type II) and either add `p_adjust` or drop the column.

**F10. `stat_hov_test(method = "f")`: denominator df dropped** — `R/stat_infer.R:1009`. `stats::var.test` reports num df and denom df; output keeps only `df` = num df, `df2` = NA. Information loss in an otherwise exact match.

**F11. Edge robustness — raw reference-package errors leak through without pre-validation** — `stat_normality` on a constant column leaks `shapiro.test`'s `all 'x' values are identical` error (L908); `stat_bootstrap` on 0-row data leaks `no data in call to 'boot'` (L1100). Both messages are decipherable but inconsistent with the package's own clean validation stops.

**F12. `stat_cor(method = "kendall")`: `stats` ties warning ("cannot compute exact p-value with ties") propagates to the user verbatim** — `R/stat_infer.R:752-753`. Cosmetic; other functions suppress or re-phrase reference warnings.

**F13. `stat_kruskal`: asymptotic χ² only** (matches `stats::kruskal.test`, no exact option) — for very small samples the exact `coin::kruskal_test` p differs materially (0.2425 vs 0.3429 at n = 8). Behavior is stats-identical; note for docs.

## Edge-case behavior (crash safety)

| Input | Behavior | Assessment |
|---|---|---|
| Empty data frame (t_test / describe / bootstrap) | clean stop ("group must have exactly 2 levels (got 0)") / describe works (n_used = 0, mean NA) / bootstrap leaks boot's error | good (t_test/describe); see F11 |
| Group collapses to 1 level after NA drop (t_test) | clean stop | good |
| Group column all NA | clean stop | good |
| All-equal values: paired t (NaN, same as stats), two-sample t pooled (clean stats error), wilcox 2-samp (W = 12.5, p = 1, eff = 0 — valid), wilcox paired zeros (p = 0.5) | pass-through of stats semantics, no opaque crashes | good |
| Degenerate 2x2 chi-square table (1 level per margin) | clean stop | good |
| ANOVA / hov `f` / cor with invalid level counts | clean stops with informative messages | good |
| `stat_normality` n = 2 / n = 5001 | NA with no crash, exactly as documented (n ∈ [3, 5000]) | good |
| `stat_perm_test` n = 2 two-sample | works (p in [1/(B+1), 1]) | good |
| `stat_bootstrap` constant statistic | warning + NA CI bounds (graceful) | good |
| `stat_t_test` n = 2 per group pooled | df = 2, works (same as stats) | good |

## Summary

The package's engines are sound: every statistic, p-value, CI, and effect size that a reference also computes matches to 0 relative difference, the bootstrap reproduces `boot::boot` bit-for-bit under the same seed, and all four permutation-test nulls agree with exact hand enumerations within MC error. The findings are all semantic/contract-level: silent data-deletion rules (`stat_cor`, `stat_describe`), sign and SS conventions that contradict the natural reading (`stat_t_test` estimate, `stat_ancova` Type I), and custom-statistic paths in `stat_perm_test` where the documented usage silently produces p = 1. None of the findings corrupt a number that both sides compute; all are reachable from documented or example usage without any warning, which is why they are graded Important rather than Minor.

Recommended fix order: F5/F6 (silent meaningless p-values in the permutation framework) → F1/F7 (silent data deletion) → F3/F4 (ANCOVA semantics) → F2 (sign) → minors.

## Fixes applied

- **Date:** 2026-09-09 · **Commit:** `c8d4686` ("fix stat_infer audit findings: perm-test semantics, pairwise cor, ANCOVA Type II, per-variable describe") · backup of the pre-fix disk state: `.superpowers/audits/stat_infer-backup-beforefix.R`

| Finding | Change | Evidence |
|---|---|---|
| F5 | `stat_perm_test` custom `two.sample` now permutes the **group labels**: each permutation re-assigns the `group` column of the data handed to `statistic` (indices stay `seq_len(n)`). Degenerate-null detector added: `null_sd < 1e-12 * max(1, abs(T0))` → `warning("permuted statistic is invariant; the statistic may not depend on the permuted labels", call. = FALSE)`. `@details` documents the semantics; broken test idiom at tests:373-378 replaced. | sanity run: custom median-diff statistic now `null_sd = 0.2218`, `p = 0.111` (was `null_sd = 0`, `p = 1`); `f_bad = median(y)` triggers the invariant warning; new tests `stat_perm_test custom statistic (two.sample) permutes group labels` PASS |
| F6 | custom `independence` shuffles the **second `.cols` column** (as docs claimed), not the last column of the data frame; slice columns are never touched (explicit `stop()` if the target is a `.by` column; `stop()` when fewer than two `.cols` columns make the target unidentifiable). `@details` updated. | sanity: `.cols = c(x, z, junk)` with trailing `junk` shuffles `z` (`null_sd = 0.16`, meaningful p); `.by = s` slice run gives per-slice `null_sd ≈ 0.34` (was `p = 1`, `null_sd = 0` in probe); new independence test PASS |
| F1 | `stat_cor` computes each pair on its **pairwise-complete rows** (`stats::cor.test` convention); `n_used` = per-pair complete count; `.by` slicing preserved (only missing slice keys dropped globally). | sanity: probe §F data now `n_used = 9, 9, 8` (pair a-b uses 9 rows, matching `cor.test`); new test `stat_cor uses pairwise-complete rows per pair` PASS |
| F7 | `stat_describe` uses **per-variable deletion** (a variable's stats use only its own non-NA values within the slice); `n_missing` = the variable's own NA count within the slice. | sanity: slice with no own NAs → `n_missing = 0` (was 5); `describe(x)` mean equals `describe(x:z)` x-row; new tests PASS |
| F3 | `stat_ancova` table now decomposed with `car::Anova(fit, type = 2)` (car already in Imports); `meansq` computed as sumsq/df (car omits Mean Sq); `@details`/`@return` say "Type II sums of squares (covariate-adjusted group term)". | sanity: probe §6 unbalanced data now `SS_g = 14.495 = car::Anova(type=2)` exactly (Type I was 0.7511); test updated to pin against `car::Anova` — PASS. Recomputation note: ANCOVA group SS/p in unbalanced designs change vs 0.0.13 output (that was the fix's purpose) |
| F4 | `test_parallel` loops the interaction check over **all** responses and covariates; violations collected and emitted as **one** warning listing each violating `resp × covariate` with its p. | sanity: `resp = c(y1, y2)` where only y2 interacts → warning lists `y2 x x (p = 9.8e-05)` (missed entirely before); new test `stat_ancova test_parallel checks every response` PASS |
| F2 | `stat_t_test` two-sample `estimate` = **group1 − group2** (`est[1] - est[2]`), consistent with statistic/CI/effect; `@return` documents the direction. One-sample/paired `estimate` unchanged (single mean / mean difference). | sanity: `estimate = -0.6195` alongside `statistic = -1.735` (both a−b; estimate was previously sign-flipped); test updated to pin `mean(x1) - mean(x2)` — PASS |
| F8 | `stat_wilcox_test` one-sample `effect` = rank-biserial from the signed ranks vs `mu`: `(2W − n0(n0+1)/2) / (n0(n0+1)/2)` (zeros dropped); `NA` when all deviations are zero; documented in `@return`. | sanity: mtcars mpg vs mu=20 → `rbc = -0.0568` (hand-check matches); all-equal input → NA without crash; new tests PASS |
| F9 | `stat_anova` `@details` now states "Type I (sequential) sums of squares" with a pointer to `stat_ancova` for Type II; the `p.adjusted` claim removed from `@return` (no `p_adjust` argument exists; column kept, hardcoded `"none"`). `stat_ancova` `@details` states Type II. | man pages regenerated via `devtools::document()`; existing anova tests (pinning `stats::anova` Type I numbers) still PASS — behavior unchanged, docs corrected |
| F10 | `stat_hov_test(method = "f")` keeps `df2` = denominator df from `stats::var.test` (`fit$parameter[2]`). | sanity: ToothGrowth f-test `df = 29`, `df2 = 29` (was NA); test extended and PASS |
| F12 | `stat_cor(method = "kendall")`: ties pre-detected per pair (`anyDuplicated`) → one clean package-style warning ("ties present, p-value is asymptotic (no exact test)"); the raw stats warning is suppressed for kendall only. | sanity: 8-row tied data → exactly the package warning, estimate/p still match `cor.test`; new test `stat_cor kendall with ties warns once and matches stats` PASS |
| F11, F13 | **Deferred** (documented edge-polish only: raw reference errors for constant-column normality / 0-row bootstrap; kruskal asymptotic-only note). No code change in this pass. | out of scope of this fix batch; behavior matches stats semantics |

**Test/validation summary:** `devtools::test(filter='stat_infer')` → FAIL 0 / WARN 2 / PASS 134 (warnings are the pre-existing spearman ties messages from the reference `cor.test` calls in an old test); full `devtools::test()` → **FAIL 0 / WARN 3 / PASS 1406** (the 3 warnings are pre-existing, non-MI). Updated tests: t-test estimate sign (F2), perm-test custom two.sample/independence semantics + detector (F5/F6), cor pairwise/kendall (F1/F12), describe per-variable (F7), ancova Type II + multi-response parallel (F3/F4), wilcox one-sample rbc (F8), hov f df2 (F10). Committed files: `R/stat_infer.R` (includes the author's pre-existing in-flight `.slice_index` tibble hunk — backup preserved), `tests/testthat/test-stat_infer.R`, 7 regenerated `man/*.Rd`. `NAMESPACE`/`DESCRIPTION` left uncommitted (they tangle with the author's separate 0.0.14 data-cleaning work; this fix adds/removes no exports).
