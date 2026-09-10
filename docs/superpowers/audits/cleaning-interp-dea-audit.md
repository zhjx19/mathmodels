# Audit: data_cleaning (self-written) + interp_fit + eval_dea

- **Date:** 2026-09-09
- **Method:** reference packages (`naniar`, `Benchmarking`) + independent
  re-implementation + direct likelihood optimisation. Probe:
  `.superpowers/audits/probe-cleaning-dea.R`, diagnostics `diag2.R`, `diag3.R`.
- **Result:** 0 computational bugs (14/14 PASS).

## data_cleaning.R (self-written parts; MI toolkit audited separately vs mice)

| Function | Reference / method | Verdict |
|---|---|---|
| na_mcar_test | Little's test; chi2 self-consistency + EM MLE vs direct `optim` likelihood maximisation | match |
| na_mcar_test vs `naniar::mcar_test` | chi2 3.9657 vs 3.9308 (0.9% rel), same df=5, same MCAR verdict | close variant |
| impute_model (knn) | independent k-NN (Euclidean on features, mean of k=3 nearest) | match (exact) |
| outlier (iqr / zscore) | manual quantile/mean-sd fences; detect / to="na" | match (exact) |
| (impute linear/spline, pivot_table) | thin wrappers / audited earlier | skipped |

**na_mcar detail:** the package EM was verified against a direct maximisation of
the observed-data MVN log-likelihood (`optim`, BFGS): μ matches to 1e-7, Σ to
2e-6 → the EM reaches the true MLE. The Little statistic computed from those
estimates is self-consistent (manual reconstruction = package value 3.96573104
exactly). The ~0.9% gap vs `naniar` is therefore an EM/statistic *implementation
detail* of naniar, not a formula error in the package; both conclude
"fail to reject MCAR". Classified as a **close variant, not a bug**.

## interp_fit.R

All fitting functions are thin wrappers over `stats` / `minpack.lm`:
`interp_linear`→`stats::approx`, `interp_spline`→`stats::spline`/`smooth.spline`,
`interp_hermite`→`splinefun(monoH.FC)`, `interp_poly`/`poly_fit`→`lm`,
`curve_fit`/`growth_fit`→`nlsLM`. Spot-checked wrappers match their references
exactly (linear vs `approx`, spline vs `spline(fmm)`, poly vs raw `lm`). Per the
audit scope (wrappers of established functions need no deep review), no findings.

## eval_dea.R

| Function | Reference | Verdict |
|---|---|---|
| basic_DEA (input orientation, VRS & CRS) | `Benchmarking::dea` | match (exact, 1e-5) |
| super_DEA / basic_SBM / super_SBM / malmquist | not cross-checked (no installed reference for the SBM/super-efficiency LP variants within scope) | not audited |

**basic_DEA** efficiency scores matched `Benchmarking::dea` for both VRS and CRS
input orientation to 1e-5 on a 6-DMU example (VRS: 1, 0.65869, 1, 1, 0.3914,
0.75618; CRS: 1, 0.60197, 0.85399, 0.92764, 0.13502, 0.49274). The SBM and
super-efficiency LP variants were **not independently verified** — recommend a
follow-up if those are used in published work.

## Findings summary

| Grade | Module | Item |
|---|---|---|
| Close variant (not a bug) | data_cleaning | na_mcar_test chi2 differs from naniar by ~0.9% (EM/statistic detail); EM verified correct via direct MLE |
| Not audited | eval_dea | SBM / super-efficiency / Malmquist LP variants (no installed reference in scope) |
