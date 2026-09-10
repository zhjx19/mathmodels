# Audit: pred_grey + eval_inequality (hand-written, no/partial reference packages)

- **Date:** 2026-09-09
- **Method:** independent re-implementation of each textbook formula (coded from the
  literature, NOT copied from the package) + reference packages where they exist
  (`ineq`, `DescTools`, `pracma`). The CRAN `grey` package is withdrawn, so
  pred_grey is closed-form math only.
- **Probes:** `.superpowers/audits/probe-pred_grey.R`, `.superpowers/audits/probe-inequality.R`
- **Result:** pred_grey **31/31 PASS, 0 bugs**; eval_inequality **1 real bug found & fixed**
  (commit `f00454f`), 21/21 PASS after the fix.

## pred_grey.R — GM11, GM1N, DGM21, verhulst, combine_preds

| Function | Check | Verdict |
|---|---|---|
| GM11 | params a/b vs independent OLS; fitted/pnext vs independent time response; f(n+1)=pnext; exact-ish recovery on `c·e^{λk}` (rel err 3.9e-4); constant/n<4/negative all stop cleanly | match |
| DGM21 | params vs independent OLS; fitted vs independent closed-form reconstruction (1.6e-14); pnext; discrete-relation theory (drop boundary pair recovers `a=e^{a_c}-1`); fit quality | match |
| verhulst | params vs independent OLS; fitted/pnext vs independent closed form (1.8e-14); near-true param recovery (a=0.4904 vs 0.50) | match |
| GM1N | fitted vs independent closed form (2.3e-8, names-attribute-only difference); C and P posterior-variance metrics; new_data prediction path | match |
| combine_preds | principal eigenvector of the cosine-similarity matrix; weighted mean; n=2 and all-identical equal weights; non-negative weights | match |

**Notes (not bugs):** grey models (GM11/DGM21/Verhulst) exhibit the well-known
discretization gap between the difference-equation OLS parameters and the continuous
time-response; exact recovery is therefore not expected (GM(1,1) rel err ~4e-4;
DGM(2,1) larger). The package is a faithful implementation of the textbook method.

## eval_inequality.R

| Function | Reference / method | Verdict |
|---|---|---|
| trapz | `pracma::trapz` | match (exact) |
| gini0 | `ineq::Gini` | match (exact) |
| gini0 vs `DescTools::Gini` | convention: DescTools defaults to the n/(n-1) sample-corrected estimator; `gini0` uses the population formula | variant (documented) |
| gini (grouped) | independent Lorenz-curve trapezoid; singleton-group == gini0 | match |
| theil0 | `ineq::Theil` (Theil T) | match (exact) |
| theil (grouped) | independent implementation; singleton == theil0 | match |
| theil0_g / theil_g | identity `T == Tb + Tw`; `theil0_g` T == `theil0(micro)`; `theil_g` T == `theil(grouped)` | match |
| theil_g2_nest | identity `T == Tb_g1 + Tb_g2 + Tw`; T == `theil(raw rows)` | match |
| **theil_g2_cross** | identity `T == Tb + Tw_b + Tw_w` held, but **T != direct Theil of the data** | **BUG (fixed)** |

### Bug: `theil_g2_cross` within-group2 weight used population share, not income share

`R/eval_inequality.R` (~line 297): the within-group2 Theil of each cross cell was
weighted by the cell's **population share within group1** (`sum(dr$pop)/sum(d$pop)`)
instead of its **income share within group1** (`sum(dr$y*dr$pop)/sum(d$y*d$pop)`).
Theil decomposition weights are always income shares, so the cross decomposition did
not add up to the true total.

Reproduction (doc example `data3`):
```
package TT                    = 0.27983121
theil(raw rows) (true total)  = 0.27882417     # mismatch
independent sequential decomp = 0.27882417     # exact reconstruction
```
The nested decomposition (`theil_g2_nest`) correctly used income shares and was
unaffected; `theil0_g`/`theil_g` were also correct.

**Fix** (commit `f00454f`): weight each cell's within-Theil by its income share
within group1. Added regression tests asserting `theil_g2_cross` (and
`theil_g2_nest`) satisfy `T == direct Theil` and the component-sum identity.
After the fix: `test-inequality` 16/16 pass; probe 21/21 pass.

## Findings summary

| Grade | Count | Items |
|---|---|---|
| Critical | 1 | `theil_g2_cross` wrong within-component weight → decomposition total wrong (fixed) |
| Important | 0 | — |
| Minor | 0 | (grey discretization bias + DescTools Gini convention are method/convention notes, not defects) |
