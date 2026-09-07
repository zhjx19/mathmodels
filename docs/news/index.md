# Changelog

## mathmodels 0.0.12

### Breaking changes

- **Renamed** the differential equation model functions from `model_*()`
  to `ode_*()` for a consistent prefix with
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md):
  `model_malthus()` →
  [`ode_malthus()`](https://zhjx19.github.io/mathmodels/reference/ode_malthus.md),
  `model_logistic()` →
  [`ode_logistic()`](https://zhjx19.github.io/mathmodels/reference/ode_logistic.md),
  `model_si()` →
  [`ode_si()`](https://zhjx19.github.io/mathmodels/reference/ode_si.md),
  `model_sis()` →
  [`ode_sis()`](https://zhjx19.github.io/mathmodels/reference/ode_sis.md),
  `model_sir()` →
  [`ode_sir()`](https://zhjx19.github.io/mathmodels/reference/ode_sir.md),
  `model_seir()` →
  [`ode_seir()`](https://zhjx19.github.io/mathmodels/reference/ode_seir.md),
  `model_lv()` →
  [`ode_lv()`](https://zhjx19.github.io/mathmodels/reference/ode_lv.md).
- **Removed** the redundant `method` argument of
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md)
  (the two methods were numerically identical);
  `R_t = beta * S / (gamma * N)` is now always used and `N` is inferred
  from the known compartments (S, E, I, R) only.

### Bug fixes

- **pred_reg**: fixed the Hosmer-Lemeshow statistic (wrong column
  reference and missing variance term — previously reported absurd
  chi-squares); fixed the factor-response coding in
  [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md)
  (events were mapped to -1/0); logistic prediction CIs now use the link
  scale (previous CIs were far too narrow); `.predict_lm()` returns both
  `confidence_interval` and `prediction_interval` (the prediction
  interval was mislabeled);
  [`reg_predict()`](https://zhjx19.github.io/mathmodels/reference/reg_predict.md)
  accepts a `newdata` argument and returns the data used; rank-deficient
  models no longer break the coefficient table;
  [`shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html) is
  guarded for n \> 5000.
- **eval_algorithms**:
  [`topsis()`](https://zhjx19.github.io/mathmodels/reference/topsis.md)
  validates `index` and `w` and rejects all-zero columns; sample names
  default to `Sample1..n`;
  [`rank_sum_ratio()`](https://zhjx19.github.io/mathmodels/reference/rank_sum_ratio.md)
  rejects unknown `method`;
  [`grey_corr()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  guards the degenerate all-equal case and validates `w`.
- **eval_dea**:
  [`super_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  now implements Tone (2002) super-efficiency (efficient DMUs get rho
  \>= 1 instead of NA); SBM undesirable-output weights follow
  Tone (2004) (1/s1, 1/s2); `basic_DEA(ud_outputs=)` now actually treats
  undesirable outputs as inputs (Hailu-Veeman) instead of silently
  ignoring them; `malmquist(type1 = "seq")` implements sequential
  reference sets (was identical to “cont”); character periods sort
  numerically; strictly positive inputs/outputs are enforced; SBM
  rejects the inapplicable `orientation` argument; malmquist warns on
  panel irregularities.
- **eval_weights**:
  [`critic_weight()`](https://zhjx19.github.io/mathmodels/reference/critic_weight.md)
  /
  [`entropy_weight()`](https://zhjx19.github.io/mathmodels/reference/entropy_weight.md)
  accept matrices again (matrix inputs previously errored element-wise);
  non-negative data is enforced for entropy paths;
  [`combine_weights()`](https://zhjx19.github.io/mathmodels/reference/combine_weights.md)
  rejects unknown `type`; `game_linear()` falls back to a pseudo-inverse
  for proportional weights;
  [`AHP()`](https://zhjx19.github.io/mathmodels/reference/AHP.md)
  returns CI = CR = 0 for n \<= 2 and rejects non-positive matrices;
  unknown `index` values are rejected instead of silently treated as
  “already normalized”.
- **eval_fuzzy**:
  [`defuzzify()`](https://zhjx19.github.io/mathmodels/reference/defuzzify.md)
  rejects unknown methods and returns the first maximum on ties;
  [`fuzzy_eval()`](https://zhjx19.github.io/mathmodels/reference/fuzzy_eval.md)
  validates R in \[0, 1\] with no NA; membership functions validate
  their parameters (strict ordering / positivity).
- **pred_grey / pred_ts**:
  [`combine_preds()`](https://zhjx19.github.io/mathmodels/reference/combine_preds.md)
  handles identical predictions (equal weights); GM models reject
  non-positive/NA data and near-constant series (a ~ 0) with informative
  errors;
  [`GM1N()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  returns `C`/`P` as documented;
  [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md)
  applies Box-Cox bias adjustment to point forecasts (`bias_adj = TRUE`
  default); Ljung-Box diagnostics use ARMA-only fitdf and no longer
  produce NaN rows for short series;
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
  /
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md)
  validate differencing depth and duplicated levels.
- **diffEq**:
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  validates `times` (strictly increasing) and rejects unused
  compartments in `init`; equation parse errors name the offending
  variable;
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md)
  warns on non-uniform time steps;
  [`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md)
  infers N from known compartments only and validates beta/gamma; the
  `b` parameter is documented as constant recruitment.
- **eval_system / eval_regional / eval_preprocess / interp_fit**:
  coupling degree validates `type` and unit-sum weights; `id_cols` typos
  are rejected;
  [`HHI()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  handles single values and rejects negatives; entropy / CV weight
  functions reject NA; preprocessing helpers guard degenerate inputs
  (`rescale_middle`, `to_positive`, `rescale_extreme`, `normalize`,
  `standardize`);
  [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md)/[`growth_fit()`](https://zhjx19.github.io/mathmodels/reference/growth_fit.md)
  report AIC/BIC on the original scale (fixes NaN at n = p) with
  non-finite start-value fallbacks;
  [`interp_poly()`](https://zhjx19.github.io/mathmodels/reference/interp_poly.md)
  validates `degree`.
- **read_nbs**: `range = "A1"` (version-robust header read).
- **DESCRIPTION**: `dplyr (>= 1.1.0)` now declared
  ([`reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
  requires it).

### Multivariate statistics toolkit (new)

- **Dimension reduction** (`mv_stats.R`):
  - [`mv_pca()`](https://zhjx19.github.io/mathmodels/reference/mv_pca.md):
    Principal component analysis with tidy scores, loadings, and
    variance decomposition.
  - [`mv_efa()`](https://zhjx19.github.io/mathmodels/reference/mv_efa.md):
    Exploratory factor analysis (maximum likelihood) with loadings,
    uniquenesses, variance shares, and optional factor scores.
- **Clustering** (`mv_stats.R`):
  - [`mv_hclust()`](https://zhjx19.github.io/mathmodels/reference/mv_hclust.md):
    Hierarchical clustering with tidy cluster assignments, merge
    heights, and the raw `hclust` object for custom plots.
  - [`mv_kmeans()`](https://zhjx19.github.io/mathmodels/reference/mv_kmeans.md):
    K-means clustering with assignments, centers, and within-cluster
    sums of squares.
  - [`mv_cluster_metrics()`](https://zhjx19.github.io/mathmodels/reference/mv_cluster_metrics.md):
    Clustering validity metrics computed in base R: per-observation
    silhouette widths, average silhouette, Davies-Bouldin,
    Calinski-Harabasz, and Dunn indices.
  - [`plot_mv_elbow()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_elbow.md)
    and
    [`plot_mv_silhouette()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_silhouette.md):
    Elbow plot (k-means or hierarchical) and silhouette plot for both
    clustering types;
    [`plot_mv_dendrogram()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_dendrogram.md)
    for hierarchical dendrograms.
- **Discriminant analysis** (`mv_stats.R`):
  - [`mv_discrim()`](https://zhjx19.github.io/mathmodels/reference/mv_discrim.md):
    LDA/QDA via `MASS` with training predictions, confusion matrix,
    accuracy, discriminant scores (LDA), and optional classification of
    new data.
- **Other multivariate methods** (`mv_stats.R`):
  - [`mv_corresp()`](https://zhjx19.github.io/mathmodels/reference/mv_corresp.md):
    Simple correspondence analysis with row/column coordinates and
    inertia decomposition.
  - [`mv_mds()`](https://zhjx19.github.io/mathmodels/reference/mv_mds.md):
    Metric (`cmdscale`) and non-metric
    ([`MASS::isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html))
    multidimensional scaling.
  - [`mv_cancor()`](https://zhjx19.github.io/mathmodels/reference/mv_cancor.md):
    Canonical correlation analysis with Wilks’ lambda significance
    tests, coefficients, and structure correlations.
  - [`mv_pcor()`](https://zhjx19.github.io/mathmodels/reference/mv_pcor.md):
    Partial correlation matrix with t-tests.
- **Visualization**:
  [`plot_mv_scree()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_scree.md)
  and
  [`plot_mv_scores()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_scores.md)
  for PCA/EFA scree and score plots,
  [`plot_mv_corresp()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_corresp.md)
  and
  [`plot_mv_mds()`](https://zhjx19.github.io/mathmodels/reference/plot_mv_mds.md)
  for CA/MDS configuration maps.
- **Dependencies**: unchanged. Tidy-select column interfaces are
  implemented through the existing `dplyr` dependency.
- All functions follow a tidy interface: tidy-select columns in, tibble
  components out, and [`print()`](https://rdrr.io/r/base/print.html)
  methods for compact summaries.

## mathmodels 0.0.11

- **New function**
  [`ts_arimax()`](https://zhjx19.github.io/mathmodels/reference/ts_arimax.md):
  Fits ARIMA models with exogenous regressors (ARIMAX). Supports both
  automatic order selection
  ([`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html))
  and manual specification
  ([`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html)),
  consistent with the existing
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md)
  interface.
- **Updated**
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md):
  Added `newxreg` parameter for forecasting ARIMAX models with future
  exogenous regressor values.
- **Improved**
  [`compute_mf()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md)
  /
  [`compute_mf_funs()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md):
  Added `.builder` parameter supporting Gaussian (`"gauss"`) and sigmoid
  (`"sigmoid"`) membership functions in addition to the default
  triangular/trapezoidal (`"tri"`). Custom builder functions are also
  accepted. Renamed `thresholds` parameter to `knots` for clarity. Fixed
  [`defuzzify()`](https://zhjx19.github.io/mathmodels/reference/defuzzify.md)
  floating-point comparison in `max_membership` method.

### Interpolation & curve fitting toolkit (new)

- **Interpolation** (`interp_fit.R`):
  - [`interp_linear()`](https://zhjx19.github.io/mathmodels/reference/interp_linear.md):
    Piecewise linear interpolation.
  - [`interp_spline()`](https://zhjx19.github.io/mathmodels/reference/interp_spline.md):
    Cubic spline interpolation with optional smoothing via `spar`.
  - [`interp_poly()`](https://zhjx19.github.io/mathmodels/reference/interp_poly.md):
    Polynomial interpolation of specified degree.
  - [`interp_hermite()`](https://zhjx19.github.io/mathmodels/reference/interp_hermite.md):
    Shape-preserving piecewise cubic Hermite interpolation.
- **Curve fitting** (`interp_fit.R`):
  - [`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md):
    Polynomial regression.
  - [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md):
    Unified interface for linearizable curves via variable
    transformation (exponential, power-law, logarithmic, hyperbolic).
  - [`growth_fit()`](https://zhjx19.github.io/mathmodels/reference/growth_fit.md):
    Nonlinear growth curve fitting via
    [`minpack.lm::nlsLM()`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html)
    with automatic starting values (Logistic, Gompertz, exponential
    saturation, Michaelis-Menten).
- **Refinements**:
  - [`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md)
    now returns explicit polynomial equations (e.g.,
    `y = 3 + 2 * x - 0.5 * x^2`) instead of R formula strings.
  - Internal `._compute_fit_stats()` avoids redundant AIC/BIC
    computation for `lm` models (uses
    [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) /
    [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) directly;
    normal-approximation fallback for `nls`).
  - [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md)
    documentation now notes log-transform back-transformation bias for
    `"exp"` and `"power"` types.
- **Dependencies** (new): `minpack.lm` for Levenberg-Marquardt nonlinear
  least squares.

## mathmodels 0.0.10

### Regression prediction toolkit (new)

- **Model fitting** (`pred_reg.R`):
  - [`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md):
    OLS multiple linear regression with optional stepwise selection
    (forward / backward / both), returning coefficients with CIs and fit
    statistics.
  - [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md):
    Binary logistic regression (`glm(binomial)`) with stepwise
    selection, odds ratios, and Hosmer-Lemeshow goodness-of-fit test.
  - [`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md):
    Poisson regression for count data, with dispersion check.
  - [`reg_negbin()`](https://zhjx19.github.io/mathmodels/reference/reg_negbin.md):
    Negative binomial regression via
    [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) for
    over-dispersed counts.
- **Diagnostics and prediction** (`pred_reg.R`):
  - [`reg_diagnostics()`](https://zhjx19.github.io/mathmodels/reference/reg_diagnostics.md):
    Unified diagnostics — Breusch-Pagan, Durbin-Watson, Shapiro-Wilk
    (LM); Hosmer-Lemeshow, dispersion/deviance (GLM).
  - [`reg_predict()`](https://zhjx19.github.io/mathmodels/reference/reg_predict.md):
    Draw new predictor data from the empirical distribution and generate
    fitted or predicted values with 95% confidence intervals.
- **Visualization** (`pred_reg.R`):
  - [`plot_reg_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_residuals.md):
    Residual diagnostic plots (residuals vs. fitted, normal Q-Q,
    density).
  - [`plot_reg_predict()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_predict.md):
    Prediction plots with confidence bands — supports both fitted
    (training data) and predicted (new data) modes.

## mathmodels 0.0.9

### Time series toolkit (new)

- **Data structure and preprocessing** (`pred_ts.R`):
  - [`ts_df()`](https://zhjx19.github.io/mathmodels/reference/ts_df.md)
    and
    [`as_ts_df()`](https://zhjx19.github.io/mathmodels/reference/as_ts_df.md):
    Lightweight single-series time-series data frame and conversion
    helpers.
  - [`complete_ts_df()`](https://zhjx19.github.io/mathmodels/reference/complete_ts_df.md)
    and
    [`impute_ts_df()`](https://zhjx19.github.io/mathmodels/reference/impute_ts_df.md):
    Gap completion and linear / spline interpolation.
- **Modeling and diagnostics** (`pred_ts.R`):
  - [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md):
    Stationarity diagnostics (ADF and KPSS).
  - [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md):
    Seasonal-trend decomposition (STL).
  - [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md):
    Exponential smoothing.
  - [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md):
    SARIMA with auto-selection or user-specified orders.
  - [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md):
    Unified forecasting for ETS and SARIMA models.
  - [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
    and
    [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md):
    Box-Cox / log / differencing workflows and inverse transformation
    for forecasts.
- **Visualization** (`pred_ts.R`):
  - [`plot_ts()`](https://zhjx19.github.io/mathmodels/reference/plot_ts.md),
    [`plot_ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_forecast.md),
    [`plot_ts_acf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_acf.md),
    [`plot_ts_pacf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_pacf.md),
    [`plot_ts_stl()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_stl.md),
    and
    [`plot_ts_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_residuals.md).
- **Dependencies** (new): forecast (ETS, ARIMA), tseries (ADF/KPSS
  tests), patchwork (multi-panel plots).

### Epidemic visualization & metrics overhaul

- **Renamed** `epidemic_plot.R` → `epi_plots.R`. Consolidated from 3
  files into 1.
- **Removed** `compute_incidence()`, `plot_infectious_curve()`,
  `plot_cumulative_infection()`.
- **Merged** `epidemic_metrics()` →
  [`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md):
  simplified to 4 core scalars (R0, peak_infection, peak_time,
  attack_rate).
- **Refined**
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md):
  now shows ΔI (daily new infections) with peak time annotation.

### Testing

- New: `test-ts.R` (52 tests), `test-epi_plots.R` (19 tests).
- Full suite: **539 tests, 0 failures**.

## mathmodels 0.0.8

### New features

- **Differential equation models** (`diffEq_models.R`):
  - [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md):
    General-purpose string-formula ODE solver powered by **deSolve**.
  - `model_malthus()`, `model_logistic()`, `model_si()`, `model_sis()`,
    `model_sir()`, `model_seir()`, `model_lv()` — all using a unified
    `init` + `params` interface with physics-meaningful named parameters
    (e.g., `beta`, `gamma`, `sigma`).
- **Epidemic visualization** (`epidemic_plot.R`):
  - [`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md):
    Faceted/overlaid line plot of selected compartments.
  - `compute_incidence()`: Daily incidence from S or I differences.
  - [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
    `plot_infectious_curve()`, `plot_cumulative_infection()`: Infection
    curve visualizations.
  - [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md):
    S–I phase portrait.
  - [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md):
    Effective reproduction number trajectory.
- **Epidemic metrics** (`epidemic_metrics.R`):
  - `epidemic_metrics()`: Comprehensive epidemic summary (R0, peak,
    attack rate, trajectory with Rt and growth rate).
- **Markov chain prediction** (`markov.R`):
  - [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md):
    Transition probability matrix, multi-step prediction, stationary
    distribution.
  - [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md):
    Grey–Markov combined prediction with GM(1,1) and Markov correction.

### Dependencies

- **Imports** (new): Added `deSolve` for ODE solving.

### Testing

- 68 tests for diffEq, 23 for epidemic plots, 21 for epidemic metrics,
  14 for Markov chain models. Zero errors.

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
