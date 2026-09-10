# Audit: eval_weights + eval_fuzzy (hand-written formulas, no reference packages)

- **Date:** 2026-09-09
- **Method:** independent re-implementation of each textbook formula + reference
  packages where they exist (`stats::prcomp`, `stats::factanal`). Probes:
  `.superpowers/audits/probe-weights.R`, `.superpowers/audits/probe-fuzzy.R`.
- **Result:** eval_weights **1 real bug found & fixed** (AHP RI table, commit
  `ead97af`); eval_fuzzy **0 bugs** (38/38 PASS).

## eval_weights.R

| Function | Check | Verdict |
|---|---|---|
| AHP | w vs independent principal eigenvector; Lmax; CI; CR = CI/RI[n]; perfect-consistency recovers w with CR=0 | match, **except RI[12]** |
| combine_weights | linear / multiplicative / game (geometric mean) / game_linear (normal-equation solve) all vs independent formulas | match |
| weight_cv | cv = sd/mean normalized | match |
| weight_entropy | rescale to [ε,1-ε], p=x/Σx, e=−Σp·ln p/ln n, w=(1−e)/Σ(1−e) | match |
| weight_critic | std × Σ(1−r) (std and entropy contrast) | match |
| weight_pca | λ vs `prcomp` sdev²; weights sum to 1; standardization only for non-NA index (documented) | match (by design) |
| weight_fa | share = colSums(A²)/p; weights sum to 1 | match |
| linear_sum | X %*% w | match |

### Bug: AHP random consistency index for n=12

`R/eval_weights.R` (RI table): the n=12 entry was **1.48**; the canonical Saaty
RI is **1.54**. This inflates CR = CI/RI[12] (e.g. CI=0.7117 → CR=0.4809 instead
of 0.4622), so a 12-criterion matrix that is actually consistent is more likely
to be wrongly flagged inconsistent. Repro (n=12 reciprocal matrix): implied
RI = 1.48 vs canonical 1.54. **Fixed** (commit `ead97af`), with a regression test
pinning CR = CI/1.54 for n=12. Only n=12 was affected; the rest of the table
(n=1..11, 13..15) matches Saaty.

**Method/convention notes (not defects):** `weight_pca`/`weight_fa` treat
`index = NA` as "already standardized" — with the default `index = NULL` no
standardization is applied (documented, but a usability trap for raw data).
`combine_weights(type="game")` is the geometric mean and `type="game_linear"` is
the game-theoretic normal-equation method (naming documented).

## eval_fuzzy.R

| Function | Check | Verdict |
|---|---|---|
| mf_tri / mf_trap | peak/plateau values, linear pieces, one-sided ±Inf shoulders, [0,1] range | match |
| mf_gauss / mf_gbell / mf_gauss2 | formula identities, peak at centre, two-sided plateau | match |
| sigmoid / mf_sigmoid / mf_dsigmoid / mf_psigmoid | 0.5 at inflection, difference/product identities | match |
| mf_z / mf_pi / mf_s | Z/PI/S textbook shape, midpoint 0.5, range | match |
| build_mf_tri / gauss / sigmoid | shoulder + middle-level construction; compute_mf evaluation | match |
| fuzzy_eval (types 1–5) | each operator vs independent hand formula; B normalized to 1 | match |
| defuzzify | weighted_average, centroid (normalizes unnormalized μ), max_membership (first max) | match |

## Findings summary

| Grade | Module | Item |
|---|---|---|
| Important | eval_weights | AHP RI[12] = 1.48 should be Saaty 1.54 → CR wrong for n=12 (fixed `ead97af`) |
| — | eval_fuzzy | none |
