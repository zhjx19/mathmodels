# Changelog

## mathmodels 0.0.8

### New features

- **Differential equation models** (`diffEq_models.R`):
  - [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md):
    General-purpose string-formula ODE solver powered by **deSolve**.
  - [`model_malthus()`](https://zhjx19.github.io/mathmodels/reference/model_malthus.md):
    Malthusian exponential growth.
  - [`model_logistic()`](https://zhjx19.github.io/mathmodels/reference/model_logistic.md):
    Logistic population growth.
  - [`model_si()`](https://zhjx19.github.io/mathmodels/reference/model_si.md):
    SI epidemic model (Susceptible–Infected).
  - [`model_sis()`](https://zhjx19.github.io/mathmodels/reference/model_sis.md):
    SIS epidemic model.
  - [`model_sir()`](https://zhjx19.github.io/mathmodels/reference/model_sir.md):
    SIR epidemic model.
  - [`model_seir()`](https://zhjx19.github.io/mathmodels/reference/model_seir.md):
    SEIR epidemic model.
  - [`model_lv()`](https://zhjx19.github.io/mathmodels/reference/model_lv.md):
    Lotka–Volterra predator–prey model. All model functions use a
    unified `init` + `params` interface with physics-meaningful named
    parameters (e.g., `beta`, `gamma`, `sigma`).
- **Epidemic visualization** (`epidemic_plot.R`):
  - [`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md):
    Faceted or overlaid line plot of selected compartments.
  - [`compute_incidence()`](https://zhjx19.github.io/mathmodels/reference/compute_incidence.md):
    Compute daily incidence from S or I differences.
  - [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md):
    Daily new infection curve.
  - [`plot_infectious_curve()`](https://zhjx19.github.io/mathmodels/reference/plot_infectious_curve.md):
    Infectious population over time.
  - [`plot_cumulative_infection()`](https://zhjx19.github.io/mathmodels/reference/plot_cumulative_infection.md):
    Cumulative infection curve.
  - [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md):
    S–I phase portrait.
  - [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md):
    Effective reproduction number (Rt) trajectory.
- **Epidemic metrics** (`epidemic_metrics.R`):
  - [`epidemic_metrics()`](https://zhjx19.github.io/mathmodels/reference/epidemic_metrics.md):
    Comprehensive epidemic summary returning R0, peak infection/time,
    attack rate, control time steps, time above threshold, and an
    augmented trajectory data frame with Rt and growth rate columns.
- **Markov chain prediction** (`markov.R`):
  - [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov.md):
    Transition probability matrix construction, multi-step state
    prediction, and stationary distribution.
  - [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/markov.md):
    Grey–Markov combined prediction with GM(1,1) and Markov chain
    correction on relative error states.

### Dependencies

- **Imports** (new): Added `deSolve` for ODE solving.

### Documentation and testing

- Added full `@examples` sections for all new and renamed functions.
- 68 tests for `diffEq_solver`, 23 for `epidemic_plot`, 21 for
  `epidemic_metrics`, 14 for Markov chain models. All pass with zero
  errors.
- Regenerated all `.Rd` files with
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).

## mathmodels 0.0.7

### New features

- **DEA models**: New `DEA.R` module implementing Data Envelopment
  Analysis with:
  - [`basic_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md):
    Radial models (CCR/BCC) with input/output orientation support,
    returning both Shephard distances and Farrell efficiencies, plus
    slacks, targets, and lambda weights.
  - [`super_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md):
    Super-efficiency radial models, fully compatible with CRS/VRS
    settings.
  - [`basic_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md):
    Slacks-Based Measure (non-radial) with CRS/VRS support.
  - [`super_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md):
    Super-efficiency SBM for distinguishing efficient DMUs.
  - [`malmquist()`](https://zhjx19.github.io/mathmodels/reference/DEA.md):
    Malmquist productivity index supporting contemporaneous, sequential,
    and global references with FGNZ and Ray-Desli decompositions.
- All DEA models use **lpSolveAPI** as the sole LP solver, with zero
  external DEA package dependency.

### Dependencies

- **Imports** (new): Added `lpSolveAPI` for linear programming in DEA
  models.
- **Suggests** (removed): Removed `deaR` from Suggests; all DEA
  functionality is now provided natively.

### Parameter validation

- Added input validation for all DEA functions, including checks for NA
  values, correct column types, and valid orientation/rts/type
  parameters.

### Testing

- Added **30 new test cases** in `test-DEA.R`, covering:
  - Basic DEA model validation (CCR, BCC, SBM)
  - Super-efficiency models
  - Malmquist index with all 6 type1 × type2 combinations
  - Undesirable outputs handling
  - Input validation (NA detection, minimum column requirements)
- All 208 tests pass with zero errors.

### Documentation

- Added `@examples` sections to all five DEA functions with
  self-contained inline data.
- Regenerated `man/DEA.Rd` with
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).

## mathmodels 0.0.6

### Bug fixes

- **[`AHP()`](https://zhjx19.github.io/mathmodels/reference/AHP.md)**:
  Extended the Random Index (RI) table from n = 11 to n = 15, and added
  validation for n \> 15. Previously,
  [`AHP()`](https://zhjx19.github.io/mathmodels/reference/AHP.md) would
  silently return `NA` for consistency ratios with more than 11
  criteria.
- **[`grey_corr()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)**:
  Fixed incorrect default weight direction in the internal validation;
  weights are now correctly validated against `nrow(cmp)` (number of
  evaluation objects) rather than `ncol(cmp)` (number of indicators).
- **[`DGM21()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)**:
  Fixed the `list(pred = fitted, ...)` return value where `fitted` was
  undefined; corrected to `list(fitted = pred, ...)`. Also fixed the
  initialization value `x0` from a hard-coded constant to the first
  element of the input series `X[1]`.
- **[`grey_corr_topsis()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)**:
  Fixed an internal transpose error in the call to
  [`grey_corr()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  that caused dimension mismatch in non-square matrices.
- **[`z_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md),
  [`pi_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md),
  [`s_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)**:
  Rewrote the two-stage spline algorithm for all three membership
  functions. The original implementation used a single-stage
  approximation with logically unreachable second-phase conditions,
  causing values to fall outside the valid \[0, 1\] range for certain
  parameter configurations. The fix ensures proper two-stage spline
  behavior with exact midpoint values of 0.5.

### Dependencies

- **Imports** (new): Added `MASS`, `purrr`, `readxl`, `rlang`, `stats`,
  `stringr`, and `tibble` to `DESCRIPTION/Imports`. These were already
  used in code via `@importFrom` in `NAMESPACE` but were missing from
  the formal dependency declaration.
- **Suggests** (new): Added `deaR` to `DESCRIPTION/Suggests` for DEA
  model functionality.

### Parameter validation

- Added input type and dimension validation
  ([`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) checks) for
  **all exported functions**, including checks for:
  - Data frame / matrix type verification
  - Vector length consistency
  - Non-negativity / positivity of weights
  - Matrix dimension bounds
  - Parameter range constraints

### Testing

- Added **178 test cases** across **18 test files**, achieving full
  coverage of all core modules:

| Module               | Test file                           | Tests |
|----------------------|-------------------------------------|-------|
| AHP                  | `test-AHP.R`                        | 8     |
| Combine predictions  | `test-combine_preds.R`              | 5     |
| Combine weights      | `test-combine_weights.R`            | 7     |
| CRITIC weighting     | `test-critic.R`                     | 5     |
| CV weighting         | `test-cv.R`                         | 4     |
| Entropy weighting    | `test-entropy.R`                    | 7     |
| Fuzzy evaluation     | `test-fuzzy.R`, `test-fuzzy_more.R` | 25    |
| Grey analysis        | `test-grey.R`                       | 10    |
| Grey models          | `test-grey_models.R`                | 16    |
| Inequality measures  | `test-inequality.R`                 | 12    |
| Linear sum           | `test-linear_sum.R`                 | 4     |
| Membership functions | `test-membership.R`                 | 19    |
| PCA weighting        | `test-pca.R`                        | 4     |
| Preprocessing        | `test-preprocess.R`                 | 21    |
| Rank sum ratio       | `test-rsr.R`                        | 5     |
| System evaluation    | `test-system_evaluation.R`          | 9     |
| TOPSIS               | `test-topsis.R`                     | 5     |

### Documentation

- Updated function examples and documentation to reflect all parameter
  validation rules.
- Added **pkgdown website** configuration (`_pkgdown.yml`). Once pushed
  to GitHub, the reference site will be automatically built and deployed
  to <https://zhjx19.github.io/mathmodels/> via GitHub Actions.
- Fixed a typo in the GitHub installation URL within `README.md`.
- Updated `README.md` version badge from 0.0.5 to 0.0.6.
- The comprehensive online manual is available at:
  <https://zhjx19.github.io/mathmodels-book/>

## mathmodels 0.0.0.9000

- Initial development version.
