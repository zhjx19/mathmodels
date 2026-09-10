# Fix report: stat_infer audit findings (F1–F13)

- **Date:** 2026-09-09
- **Repo:** E:\Github\mathmodels
- **Commit:** `c8d4686` — "fix stat_infer audit findings: perm-test semantics, pairwise cor, ANCOVA Type II, per-variable describe"
- **Files committed:** `R/stat_infer.R`, `tests/testthat/test-stat_infer.R`, `man/stat_{ancova,anova,cor,describe,perm_test,t_test,wilcox_test}.Rd`
- **Pre-fix disk backup:** `.superpowers/audits/stat_infer-backup-beforefix.R` (preserves the author's in-flight `.slice_index` hunk state)
- **Audit doc updated:** `docs/superpowers/audits/stat_infer-audit.md` — "## Fixes applied" section appended

## Status: COMPLETE

All requested fixes (F1–F10, F12) implemented and validated. F11 and F13 deferred (documented edge-polish; see Concerns).

## Fix-by-fix detail

### F5 — stat_perm_test custom two.sample: permute GROUP LABELS (R/stat_infer.R, custom branch)
- **Problem:** reordering row indices moves `y` and `g` together → any group-contrast statistic of the documented form `f(d, i)` is invariant → `null_sd = 0`, silent `p = 1`.
- **Change:** each permutation now re-assigns the `group` column of the data handed to the statistic (`d2[[g_nm]] = d[[g_nm]][sample.int(n)]`, indices stay `seq_len(n)`), so group contrasts vary. Built-in (non-custom) path untouched. Requires a `group` column with a clear `stop()` if absent. Degenerate-null detector: `null_sd < 1e-12 * max(1, abs(T0))` → `warning("permuted statistic is invariant; the statistic may not depend on the permuted labels", call. = FALSE)` — replaces silent p = 1 (applies to custom statistics in practice; built-in nulls have analytic spread). `@details` documents the semantics. Fixed the broken idiom at tests/testthat/test-stat_infer.R:373-378 (now label-based, plus an invariant-statistic `expect_warning` case).
- **Evidence:** custom median-diff statistic: `null_sd = 0.2218`, `p = 0.111` (before: `null_sd = 0`, `p = 1`). `f_bad = \(d, i) median(d$y)` → exact detector warning. New test PASS.

### F6 — custom independence: shuffle the SECOND .cols column
- **Problem:** code shuffled the LAST data-frame column (docs said "second column"); with `.by`, the last column was the slice column → silent `p = 1`.
- **Change:** shuffle target is `names(cols_df)[2]` (the second `.cols` column), permuted within each `.by` slice (`.run_one` receives the slice's `cols_df` subset); slice columns never touched. `stop()` when `ncol(cols_df) < 2` ("independence with a custom statistic requires at least two .cols columns...") and a separate `stop()` if the shuffle target is itself a `.by` column. `@details`/`@param .cols` updated.
- **Evidence:** `.cols = c(x, z, junk)`: shuffles `z` not `junk` (`null_sd = 0.16`); `.by = s` run: per-slice `null_sd ≈ 0.34`, meaningful p (probe reported `p = 1`, `null_sd = 0`). New tests PASS (including the `expect_error` for < 2 columns).

### F1 — stat_cor: per-pair pairwise-complete rows
- **Problem:** listwise deletion across ALL `.cols` (plus `.by`) shrank every pair's data to the intersection-complete subset.
- **Change:** global row filter now drops only missing `.by` keys; each pair additionally computes on its own `complete.cases` of the two variables; `n_used` = per-pair count. Documented in `@description`/`@return`.
- **Evidence:** audit repro now yields `n_used = 9, 9, 8` (pair a-b = 9, matching `cor.test`); new test `stat_cor uses pairwise-complete rows per pair` PASS.

### F7 — stat_describe: per-variable deletion + own-NA n_missing
- **Problem:** (a) listwise across all `.cols` made x's stats depend on other columns' NAs; (b) `n_missing = nrow(data) - n` counted other-slice rows and other columns' NAs.
- **Change:** slice assignment via `.by` keys only; within a slice each variable's stats use its own non-NA values (already the per-variable `x = d[[nm]][!is.na(...)]` pattern — the fix is dropping the cross-column `complete.cases`); `n_missing = n_slice - n_used` (variable's own NA count within the slice). Documented.
- **Evidence:** audit repro: slice b (no own NAs) → `n_missing = 0` (was 5); `describe(x)` x-row mean now equals `mean(x, na.rm=TRUE)` regardless of co-selected z. New tests PASS.

### F3 — stat_ancova: Type II SS via car::Anova
- **Change:** `tab = car::Anova(fit, type = 2)` (car already in Imports). car's table lacks `Mean Sq` → `meansq = sumsq/df` computed (NA for df = 0/NA). `@description`/`@return`/`@details` say "Type II sums of squares (covariate-adjusted group term)" with a note that unbalanced designs differ from sequential Type I. **stat_anova stays Type I** — its `@details` now says "Type I (sequential) sums of squares" with a pointer to stat_ancova.
- **Recomputed numbers (report note):** ANCOVA `sumsq`/`p.value` for the group term (and covariate terms) in **unbalanced** designs change vs 0.0.13 (that is the fix — group SS is now covariate-adjusted). Balanced/orthogonal designs are numerically unchanged. Audit probe example: unbalanced `SS_g` 0.7511 (Type I) → 14.495 = `car::Anova(type=2)` exactly. Test updated from `stats::anova` pinning to `car::Anova` pinning (plus a Type-I-vs-Type-II divergence assertion); PASS.
- **Existing tests updated:** `stat_ancova matches aov with covariate` → `stat_ancova uses Type II sums of squares (car::Anova)` (reason: old test pinned Type I group SS; the finding F3 mandates the change).

### F4 — test_parallel over ALL responses
- **Change:** the interaction check loops `for (rnm in names(resp_df))` × `for (cn in cov_nms)`; all violations collected into one warning listing each `` `resp` x `cov` (p = ...) ``. Message prefix unchanged ("parallel-slopes assumption may not hold: ...").
- **Evidence:** audit repro (`resp = c(y1, y2)`, only y2 interacts) now warns naming y2 (before: silent). New test asserts the warning fires and names the violating response; PASS.

### F2 — stat_t_test estimate sign
- **Change:** `estimate = est[1] - est[2]` (group1 − group2, rstatix convention), matching statistic/CI/effect direction; `@return` documents the convention. One-sample/paired estimates unchanged (they are single means / mean differences).
- **Evidence:** same-row consistency: `estimate = -0.6195` with `statistic = -1.735` (both a−b; estimate previously +flipped). Test updated from `fit$estimate[2] - fit$estimate[1]` to `mean(x1) - mean(x2)`; PASS.

### F8 — stat_wilcox_test one-sample rank-biserial
- **Change:** `effect = (2W - n0(n0+1)/2) / (n0(n0+1)/ 2)` from the signed ranks vs `mu` (zeros dropped; same construction as the paired path); `NA` when all deviations from `mu` are zero; documented in `@return` (NA-when-undefined note).
- **Evidence:** mtcars mpg vs mu=20 → `rbc = -0.0568`, hand-check `(2W - S)/S` matches; constant-at-mu input → NA, no crash. New tests PASS.

### F9 — stat_anova docs
- **Change:** `@details` adds "Type I (sequential) sums of squares" (+ comparison note re Type II/III and pointer to stat_ancova). `p.adjusted` removed from `@return` (verified in code: no `p_adjust` argument exists and L512-era hardcode `method = "none"` makes `p.adjusted == p.value` always; column kept). stat_ancova `@details` states Type II.
- **Evidence:** man pages regenerated; anova behavior unchanged → existing anova tests PASS unchanged.

### F10 — stat_hov_test method="f" df2
- **Change:** `df2 = fit$parameter[2]` for `method = "f"` (denominator df from `stats::var.test`); levene path unchanged; bartlett stays NA.
- **Evidence:** ToothGrowth f-test → `df = 29`, `df2 = 29` (df2 was NA); test extended with df/df2/p pins; PASS.

### F12 — kendall ties warning
- **Change:** pre-detect ties per pair (`anyDuplicated` on either complete column) → one clean package-style `warning("kendall correlation for `x` / `y`: ties present, p-value is asymptotic (no exact test).", call. = FALSE)`; the raw stats warning suppressed for kendall only (pearson/spearman paths untouched — spearman's stats message remains as before, consistent with existing behavior).
- **Evidence:** 8-row tied data emits exactly the new warning once per pair; estimate/p identical to `cor.test`. New test PASS.

## Deferred findings

- **F11** (raw reference errors leak for constant-column `stat_normality` / 0-row `stat_bootstrap`): messages decipherable, behavior = pass-through of stats/boot semantics; pre-validation polish left for a later pass.
- **F13** (kruskal asymptotic χ² only): stats-identical behavior; doc note candidate only.

## Updated tests (test-stat_infer.R) and reasons

| Test | Change | Reason |
|---|---|---|
| `stat_t_test two-sample Welch matches stats::t.test` | estimate pin changed to `mean(x1) - mean(x2)`; added estimate1/estimate2 pins | F2 sign convention |
| `stat_describe` + new `... per-variable deletion and own NA counts` | added | F7 |
| `stat_wilcox_test` + new `... one-sample effect is the signed-rank rbc` | added | F8 |
| `stat_ancova matches aov with covariate` → `... uses Type II sums of squares (car::Anova)`; + unbalanced Type I≠Type II assertion | repointed reference to `car::Anova(type=2)` | F3 |
| new `stat_ancova test_parallel checks every response` | added | F4 |
| `stat_cor ...` + new `... uses pairwise-complete rows per pair`, `... kendall with ties warns once and matches stats` | added | F1 / F12 |
| `stat_hov_test ...` | added df/df2/p pins for method="f" | F10 |
| `stat_perm_test custom statistic (two.sample)` → rewritten to label-permutation semantics + invariant-detector `expect_warning` | idiom fix | F5 |
| new `stat_perm_test custom statistic (independence) shuffles the second .cols column` (incl. trailing junk column + error case) | added | F6 |

## Validation

1. `devtools::test(filter='stat_infer')` → **FAIL 0 / WARN 2 / PASS 134** (the 2 warnings are the pre-existing spearman ties messages from reference `cor.test` calls inside the old spearman test; not from mathmodels code).
2. `devtools::test()` → **FAIL 0 / WARN 3 / PASS 1406** (3 pre-existing non-MI warnings, as allowed).
3. `devtools::document()` → regenerated the 7 stat_infer man pages; no export changes.
4. Commit `c8d4686`: 9 files, +382/−84.

## Git surgery notes (user's in-flight work preserved)

- The author's uncommitted `.slice_index` tibble hunk in `R/stat_infer.R` was preserved on disk, backed up to `.superpowers/audits/stat_infer-backup-beforefix.R` BEFORE edits, and committed together with my fixes (noted in the commit message) — the whole-file commit route mandated by the task.
- `NAMESPACE` and `DESCRIPTION` were NOT committed: their diffs contain only the author's separate 0.0.14 data-cleaning/multiple-imputation work (export(impute), na_*, outlier, pivot_table, ...); my fixes add/remove no exports, so skipping them loses nothing. They remain modified in the working tree for the author's own commit.
- All other in-flight files (docs/, README, NEWS, R/data_cleaning.R, tests for the new modules, ...) untouched and unstaged.

## Concerns / documented choices

1. **Two-sided p with custom statistic** keeps the existing centering (`c0 = mean(Tstar)`) — unchanged from before; not part of any finding.
2. **F5 custom two.sample now requires `group`** (labels must exist to permute). The old broken idiom worked *without* group; the new contract errors clearly when group is missing. Documented in `@param group` and `@details`.
3. **F12 suppression scope**: `suppressWarnings` wraps only the kendall `cor.test` call; if stats ever adds a *different* meaningful warning there it would be hidden — accepted trade-off per the finding ("suppress and note in doc": the package warning replaces it).
4. **F3 recomputation impact**: any downstream user comparing 0.0.13 ANCOVA group SS/p (Type I) against 0.0.14+ (Type II) on unbalanced data will see different numbers — intended, documented, and test-pinned.
5. F11/F13 deferred (see above).
