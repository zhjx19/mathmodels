# mathmodels 0.0.11

- **New function** `ts_arimax()`: Fits ARIMA models with exogenous regressors (ARIMAX). Supports both automatic order selection (`forecast::auto.arima()`) and manual specification (`forecast::Arima()`), consistent with the existing `ts_sarima()` interface.
- **Updated** `ts_forecast()`: Added `newxreg` parameter for forecasting ARIMAX models with future exogenous regressor values.
- **Improved** `compute_mf()` / `compute_mf_funs()`: Added `.builder` parameter supporting Gaussian (`"gauss"`) and sigmoid (`"sigmoid"`) membership functions in addition to the default triangular/trapezoidal (`"tri"`). Custom builder functions are also accepted. Fixed `defuzzify()` floating-point comparison in `max_membership` method.

## Interpolation & curve fitting toolkit (new)

- **Interpolation** (`interp_fit.R`):
  - `interp_linear()`: Piecewise linear interpolation.
  - `interp_spline()`: Cubic spline interpolation with optional smoothing via `spar`.
  - `interp_poly()`: Polynomial interpolation of specified degree.
  - `interp_hermite()`: Shape-preserving piecewise cubic Hermite interpolation.

- **Curve fitting** (`interp_fit.R`):
  - `poly_fit()`: Polynomial regression.
  - `curve_fit()`: Unified interface for linearizable curves via variable transformation
    (exponential, power-law, logarithmic, hyperbolic).
  - `growth_fit()`: Nonlinear growth curve fitting via `minpack.lm::nlsLM()`
    with automatic starting values (Logistic, Gompertz, exponential saturation,
    Michaelis-Menten).

- **Refinements**:
  - `poly_fit()` now returns explicit polynomial equations (e.g., `y = 3 + 2 * x - 0.5 * x^2`) instead of R formula strings.
  - Internal `._compute_fit_stats()` avoids redundant AIC/BIC computation for `lm` models (uses `stats::AIC()` / `stats::BIC()` directly; normal-approximation fallback for `nls`).
  - `curve_fit()` documentation now notes log-transform back-transformation bias for `"exp"` and `"power"` types.

- **Dependencies** (new): `minpack.lm` for Levenberg-Marquardt nonlinear least squares.

# mathmodels 0.0.10

## Regression prediction toolkit (new)

- **Model fitting** (`pred_reg.R`):
  - `reg_lm()`: OLS multiple linear regression with optional stepwise selection
    (forward / backward / both), returning coefficients with CIs and fit statistics.
  - `reg_logistic()`: Binary logistic regression (`glm(binomial)`) with stepwise
    selection, odds ratios, and Hosmer-Lemeshow goodness-of-fit test.
  - `reg_poisson()`: Poisson regression for count data, with dispersion check.
  - `reg_negbin()`: Negative binomial regression via `MASS::glm.nb()` for
    over-dispersed counts.

- **Diagnostics and prediction** (`pred_reg.R`):
  - `reg_diagnostics()`: Unified diagnostics — Breusch-Pagan, Durbin-Watson,
    Shapiro-Wilk (LM); Hosmer-Lemeshow, dispersion/deviance (GLM).
  - `reg_predict()`: Draw new predictor data from the empirical distribution and
    generate fitted or predicted values with 95% confidence intervals.

- **Visualization** (`pred_reg.R`):
  - `plot_reg_residuals()`: Residual diagnostic plots (residuals vs. fitted,
    normal Q-Q, density).
  - `plot_reg_predict()`: Prediction plots with confidence bands — supports both
    fitted (training data) and predicted (new data) modes.

# mathmodels 0.0.9

## Time series toolkit (new)

- **Data structure and preprocessing** (`pred_ts.R`):
  - `ts_df()` and `as_ts_df()`: Lightweight single-series time-series data frame and conversion helpers.
  - `complete_ts_df()` and `impute_ts_df()`: Gap completion and linear / spline interpolation.

- **Modeling and diagnostics** (`pred_ts.R`):
  - `ts_test()`: Stationarity diagnostics (ADF and KPSS).
  - `ts_stl()`: Seasonal-trend decomposition (STL).
  - `ts_ets()`: Exponential smoothing.
  - `ts_sarima()`: SARIMA with auto-selection or user-specified orders.
  - `ts_forecast()`: Unified forecasting for ETS and SARIMA models.
  - `ts_transform()` and `ts_back_transform()`: Box-Cox / log / differencing workflows and inverse transformation for forecasts.

- **Visualization** (`pred_ts.R`):
  - `plot_ts()`, `plot_ts_forecast()`, `plot_ts_acf()`, `plot_ts_pacf()`, `plot_ts_stl()`, and `plot_ts_residuals()`.

- **Dependencies** (new): forecast (ETS, ARIMA), tseries (ADF/KPSS tests), patchwork (multi-panel plots).

## Epidemic visualization & metrics overhaul

- **Renamed** `epidemic_plot.R` → `epi_plots.R`. Consolidated from 3 files into 1.
- **Removed** `compute_incidence()`, `plot_infectious_curve()`, `plot_cumulative_infection()`.
- **Merged** `epidemic_metrics()` → `epi_metrics()`: simplified to 4 core scalars (R0, peak_infection, peak_time, attack_rate).
- **Refined** `plot_incidence()`: now shows ΔI (daily new infections) with peak time annotation.

## Testing

- New: `test-ts.R` (52 tests), `test-epi_plots.R` (19 tests).
- Full suite: **539 tests, 0 failures**.

# mathmodels 0.0.8

## New features

- **Differential equation models** (`diffEq_models.R`):
  - `ode_solver()`: General-purpose string-formula ODE solver powered by **deSolve**.
  - `model_malthus()`, `model_logistic()`, `model_si()`, `model_sis()`, `model_sir()`, `model_seir()`, `model_lv()` — all using a unified `init` + `params` interface with physics-meaningful named parameters (e.g., `beta`, `gamma`, `sigma`).

- **Epidemic visualization** (`epidemic_plot.R`):
  - `plot_compartments()`: Faceted/overlaid line plot of selected compartments.
  - `compute_incidence()`: Daily incidence from S or I differences.
  - `plot_incidence()`, `plot_infectious_curve()`, `plot_cumulative_infection()`: Infection curve visualizations.
  - `plot_phase_si()`: S–I phase portrait.
  - `plot_Rt_estimate()`: Effective reproduction number trajectory.

- **Epidemic metrics** (`epidemic_metrics.R`):
  - `epidemic_metrics()`: Comprehensive epidemic summary (R0, peak, attack rate, trajectory with Rt and growth rate).

- **Markov chain prediction** (`markov.R`):
  - `markov_chain()`: Transition probability matrix, multi-step prediction, stationary distribution.
  - `GM11_markov()`: Grey–Markov combined prediction with GM(1,1) and Markov correction.

## Dependencies

- **Imports** (new): Added `deSolve` for ODE solving.

## Testing

- 68 tests for diffEq, 23 for epidemic plots, 21 for epidemic metrics, 14 for Markov chain models. Zero errors.

# mathmodels 0.0.7

## New features

- **DEA models**: New `DEA.R` module implementing Data Envelopment Analysis with:
  - `basic_DEA()`: Radial models (CCR/BCC) with input/output orientation support, returning both Shephard distances and Farrell efficiencies, plus slacks, targets, and lambda weights.
  - `super_DEA()`: Super-efficiency radial models, fully compatible with CRS/VRS settings.
  - `basic_SBM()`: Slacks-Based Measure (non-radial) with CRS/VRS support.
  - `super_SBM()`: Super-efficiency SBM for distinguishing efficient DMUs.
  - `malmquist()`: Malmquist productivity index supporting contemporaneous, sequential, and global references with FGNZ and Ray-Desli decompositions.
- All DEA models use **lpSolveAPI** as the sole LP solver, with zero external DEA package dependency.

## Dependencies

- **Imports** (new): Added `lpSolveAPI` for linear programming in DEA models.
- **Suggests** (removed): Removed `deaR` from Suggests; all DEA functionality is now provided natively.

## Parameter validation

- Added input validation for all DEA functions, including checks for NA values, correct column types, and valid orientation/rts/type parameters.

## Testing

- Added **30 new test cases** in `test-DEA.R`, covering:
  - Basic DEA model validation (CCR, BCC, SBM)
  - Super-efficiency models
  - Malmquist index with all 6 type1 × type2 combinations
  - Undesirable outputs handling
  - Input validation (NA detection, minimum column requirements)
- All 208 tests pass with zero errors.

## Documentation

- Added `@examples` sections to all five DEA functions with self-contained inline data.
- Regenerated `man/DEA.Rd` with `devtools::document()`.

# mathmodels 0.0.6

## Bug fixes

- **`AHP()`**: Extended the Random Index (RI) table from n = 11 to n = 15, and added validation for n > 15. Previously, `AHP()` would silently return `NA` for consistency ratios with more than 11 criteria.
- **`grey_corr()`**: Fixed incorrect default weight direction in the internal validation; weights are now correctly validated against `nrow(cmp)` (number of evaluation objects) rather than `ncol(cmp)` (number of indicators).
- **`DGM21()`**: Fixed the `list(pred = fitted, ...)` return value where `fitted` was undefined; corrected to `list(fitted = pred, ...)`. Also fixed the initialization value `x0` from a hard-coded constant to the first element of the input series `X[1]`.
- **`grey_corr_topsis()`**: Fixed an internal transpose error in the call to `grey_corr()` that caused dimension mismatch in non-square matrices.
- **`z_mf()`, `pi_mf()`, `s_mf()`**: Rewrote the two-stage spline algorithm for all three membership functions. The original implementation used a single-stage approximation with logically unreachable second-phase conditions, causing values to fall outside the valid [0, 1] range for certain parameter configurations. The fix ensures proper two-stage spline behavior with exact midpoint values of 0.5.

## Dependencies

- **Imports** (new): Added `MASS`, `purrr`, `readxl`, `rlang`, `stats`, `stringr`, and `tibble` to `DESCRIPTION/Imports`. These were already used in code via `@importFrom` in `NAMESPACE` but were missing from the formal dependency declaration.
- **Suggests** (new): Added `deaR` to `DESCRIPTION/Suggests` for DEA model functionality.

## Parameter validation

- Added input type and dimension validation (`stopifnot()` checks) for **all exported functions**, including checks for:
  - Data frame / matrix type verification
  - Vector length consistency
  - Non-negativity / positivity of weights
  - Matrix dimension bounds
  - Parameter range constraints

## Testing

- Added **178 test cases** across **18 test files**, achieving full coverage of all core modules:

| Module | Test file | Tests |
|--------|-----------|-------|
| AHP | `test-AHP.R` | 8 |
| Combine predictions | `test-combine_preds.R` | 5 |
| Combine weights | `test-combine_weights.R` | 7 |
| CRITIC weighting | `test-critic.R` | 5 |
| CV weighting | `test-cv.R` | 4 |
| Entropy weighting | `test-entropy.R` | 7 |
| Fuzzy evaluation | `test-fuzzy.R`, `test-fuzzy_more.R` | 25 |
| Grey analysis | `test-grey.R` | 10 |
| Grey models | `test-grey_models.R` | 16 |
| Inequality measures | `test-inequality.R` | 12 |
| Linear sum | `test-linear_sum.R` | 4 |
| Membership functions | `test-membership.R` | 19 |
| PCA weighting | `test-pca.R` | 4 |
| Preprocessing | `test-preprocess.R` | 21 |
| Rank sum ratio | `test-rsr.R` | 5 |
| System evaluation | `test-system_evaluation.R` | 9 |
| TOPSIS | `test-topsis.R` | 5 |

## Documentation

- Updated function examples and documentation to reflect all parameter validation rules.
- Added **pkgdown website** configuration (`_pkgdown.yml`). Once pushed to GitHub, the reference site will be automatically built and deployed to <https://zhjx19.github.io/mathmodels/> via GitHub Actions.
- Fixed a typo in the GitHub installation URL within `README.md`.
- Updated `README.md` version badge from 0.0.5 to 0.0.6.
- The comprehensive online manual is available at: <https://zhjx19.github.io/mathmodels-book/>

# mathmodels 0.0.0.9000

- Initial development version.
