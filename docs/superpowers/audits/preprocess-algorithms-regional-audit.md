# Audit: preprocess + algorithms + system + regional

- **Date:** 2026-09-09
- **Method:** independent re-implementation of each textbook formula (TOPSIS, GRA,
  RSR, standard coupling, EG, HHI, LQ) + manual boundary checks. Probes:
  `.superpowers/audits/probe-preprocess.R`, `probe-algo-sys-reg.R`.
- **Result:** 0 computational bugs. eval_preprocess **20/20**; algorithms/system/
  regional **16/16**. One documentation fix (commit `26fa4fc`).

## eval_preprocess.R (20/20)

| Function | Check | Verdict |
|---|---|---|
| standardize | vs `base::scale`, center/scale toggles, zero-variance stop | match |
| normalize | L2 norm, NA preserved, zero-norm stop | match |
| rescale | min-max both directions, custom [a,b], zero-range stop | match |
| rescale_middle | optimal value → 1, farthest → 0 | match |
| rescale_interval | inside interval = 1, outside linear | match |
| rescale_extreme | code `+`: x/max, `-`: min/x (standard GRA) | match (docs fixed) |
| rescale_initial / rescale_mean | x/x₁, x/mean | match |
| to_positive | minmax = max−x, reciprocal = 1/x | match |

**Doc fix (`26fa4fc`):** `rescale_extreme`'s `@return`, `@details`, and example
comments stated the *opposite* of the code (claimed `+`: min/x, `-`: x/max). The
code was correct (standard grey-relational normalization); the docs were swapped
and are now fixed.

## eval_algorithms.R (vs independent implementations)

| Function | Check | Verdict |
|---|---|---|
| topsis | L2 normalization, weighted matrix, ± ideals by `index`, Euclidean distances, S0/(S0+S*) | match (exact) |
| rank_sum_ratio | weighted rank-sum RSR = (Σ w·rank)/(Σw·n); probit regression `lm(RSR ~ Probit)` | match |
| grey_corr | ξ = (Δmin+ρΔmax)/(Δ+ρΔmax); weighted mean over sampling points | match (exact) |
| grey_corr_topsis | rescale by direction → weighted → grey correlation to ± ideals → S*/(S0+S*) | match (structure verified) |

## eval_system.R

| Function | Check | Verdict |
|---|---|---|
| coupling_degree (standard) | C = p·(∏u)^(1/p)/Σu; T = w·u; D = √(C·T) | match (exact) |
| coupling_degree (adjusted) | all-equal → C=1; range [0,1]; correction structure C=√(C1·C2) | properties hold¹ |
| obstacle_degree | O = (1−X)·w; scaled = O/rowSum (neutral 0.5 for zero rows) | match (exact) |

¹ The "adjusted" formula follows the structure cited to Wang et al. (2021); its
exact published constants could not be independently re-derived without the
paper, so it was verified by invariants (equality case, range, monotonicity)
rather than a numeric reference. Flag for the author to confirm against the
paper if publishing.

## eval_regional.R

| Function | Check | Verdict |
|---|---|---|
| LQ | (region share of industry)/(national share); first row = national benchmark | match (exact) |
| HHI | Σ(share²); scaled = (HHI−1/n)/(1−1/n); single firm = 1 | match (exact) |
| EG | (Σ(s_r−x_r)²/(1−Σx_r²) − H)/(1−H), Ellison-Glaeser (1997) | match (exact) |

## Findings summary

| Grade | Module | Item |
|---|---|---|
| Minor | eval_preprocess | `rescale_extreme` docs reversed vs code (docs fixed `26fa4fc`) |
| Note | eval_system | "adjusted" coupling formula not independently confirmable without Wang et al. (2021) |
