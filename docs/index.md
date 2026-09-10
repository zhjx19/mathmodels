# mathmodels: A Powerful R Toolkit for Mathematical Modeling

An R package providing a versatile toolkit for mathematical modeling,
developed as a companion to the book *Mathematical Modeling: Algorithms
and Programming Implementation* (China Machine Press). It focuses on
implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.14)** introduces a data-cleaning toolkit —
missing-value diagnostics with Little’s MCAR test
([`na_mcar_test()`](https://zhjx19.github.io/mathmodels/reference/na_mcar_test.md)),
missingness visuals
([`plot_na_heatmap()`](https://zhjx19.github.io/mathmodels/reference/plot_na_heatmap.md),
[`plot_na_bar()`](https://zhjx19.github.io/mathmodels/reference/plot_na_bar.md)),
simple and model-based imputation
([`impute()`](https://zhjx19.github.io/mathmodels/reference/impute.md),
[`impute_model()`](https://zhjx19.github.io/mathmodels/reference/impute_model.md)
with kNN / linear regression / decision trees), multiple imputation by
FCS chained equations
([`impute_multiple()`](https://zhjx19.github.io/mathmodels/reference/impute_multiple.md))
with Rubin’s-rules pooling
([`rubin_pool()`](https://zhjx19.github.io/mathmodels/reference/rubin_pool.md),
[`pool_fit()`](https://zhjx19.github.io/mathmodels/reference/pool_fit.md);
Barnard-Rubin small-sample df), and outlier handling
([`outlier()`](https://zhjx19.github.io/mathmodels/reference/outlier.md))
— plus the Excel-style
[`pivot_table()`](https://zhjx19.github.io/mathmodels/reference/pivot_table.md).
Together with the previous release’s tidy statistical inference toolkit
(`stat_*()` tests, bootstrap and permutation frameworks,
[`export_table()`](https://zhjx19.github.io/mathmodels/reference/export_table.md)),
the multivariate statistics toolkit (`mv_*()`), the factor-analysis
weighting
[`weight_fa()`](https://zhjx19.github.io/mathmodels/reference/weight_fa.md),
and the renamed `weight_*` / `mf_*` families. All this builds on the
regression prediction toolkit, time-series workflow, epidemic modeling,
differential equation models, grey prediction, Markov chain models,
interpolation and curve fitting, and a rich suite of evaluation methods
(AHP, Entropy, CRITIC, PCA, TOPSIS, Fuzzy, RSR, DEA).

## Key Features

- **Data Cleaning (new)** — missing-value diagnostics
  ([`na_summary()`](https://zhjx19.github.io/mathmodels/reference/na_summary.md),
  [`na_mcar_test()`](https://zhjx19.github.io/mathmodels/reference/na_mcar_test.md)
  for Little’s MCAR test) with heatmap / bar visuals; imputation: simple
  fills (group mean / median / mode, linear & spline interpolation,
  constants), a unified model-based framework
  ([`impute_model()`](https://zhjx19.github.io/mathmodels/reference/impute_model.md)
  with kNN, linear regression, decision trees), and **multiple
  imputation**
  ([`impute_multiple()`](https://zhjx19.github.io/mathmodels/reference/impute_multiple.md):
  FCS chained equations, pmm/norm for numeric columns, logreg/polyreg
  for categorical ones) with Rubin’s-rules pooling
  ([`rubin_pool()`](https://zhjx19.github.io/mathmodels/reference/rubin_pool.md)
  for scalar estimates,
  [`pool_fit()`](https://zhjx19.github.io/mathmodels/reference/pool_fit.md)
  for lm/glm fits, incl. the Barnard-Rubin small-sample df adjustment);
  every changed cell is recorded for audit. Outlier detection and
  handling
  ([`outlier()`](https://zhjx19.github.io/mathmodels/reference/outlier.md):
  IQR / z-score / Hampel-MAD bounds, set-NA or cap). Excel-style
  [`pivot_table()`](https://zhjx19.github.io/mathmodels/reference/pivot_table.md).
- **Multivariate Statistics (new)** — PCA
  ([`mv_pca()`](https://zhjx19.github.io/mathmodels/reference/mv_pca.md)),
  EFA
  ([`mv_efa()`](https://zhjx19.github.io/mathmodels/reference/mv_efa.md)),
  hierarchical and k-means clustering with validity metrics
  ([`mv_hclust()`](https://zhjx19.github.io/mathmodels/reference/mv_hclust.md),
  [`mv_kmeans()`](https://zhjx19.github.io/mathmodels/reference/mv_kmeans.md),
  [`mv_cluster_metrics()`](https://zhjx19.github.io/mathmodels/reference/mv_cluster_metrics.md)),
  discriminant analysis
  ([`mv_discrim()`](https://zhjx19.github.io/mathmodels/reference/mv_discrim.md)),
  correspondence analysis
  ([`mv_corresp()`](https://zhjx19.github.io/mathmodels/reference/mv_corresp.md)),
  MDS
  ([`mv_mds()`](https://zhjx19.github.io/mathmodels/reference/mv_mds.md)),
  canonical and partial correlation
  ([`mv_cancor()`](https://zhjx19.github.io/mathmodels/reference/mv_cancor.md),
  [`mv_pcor()`](https://zhjx19.github.io/mathmodels/reference/mv_pcor.md)),
  all with tidy components and `plot_mv_*()` visuals.
- **Statistical Inference** — tidy hypothesis tests: t-test / Wilcoxon,
  ANOVA / Kruskal-Wallis, ANCOVA, correlation, chi-square, normality and
  variance-homogeneity tests (`stat_*()`), each returning one tidy row
  per test with effect sizes; generic **bootstrap** and **permutation**
  frameworks;
  [`export_table()`](https://zhjx19.github.io/mathmodels/reference/export_table.md)
  renders results as Word/LaTeX three-line tables.
- **Evaluation Models** — AHP, Entropy weighting, CRITIC, CV, PCA and
  factor-analysis weighting, TOPSIS, Grey Relational Analysis (GRA),
  Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE) with the
  `mf_*` family of membership functions, Data Envelopment Analysis
  (CCR/BCC/SBM, Malmquist), plus inequality measures (Gini, Theil
  Index), coupling coordination degree, and obstacle degree models.
- **Prediction Models** —
  - **Grey prediction** (GM(1,1), GM(1,N), Verhulst), Markov chain
    prediction
    ([`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md),[`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)),
  - **regression prediction** toolkit
    ([`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
    [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md),
    [`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
    [`reg_negbin()`](https://zhjx19.github.io/mathmodels/reference/reg_negbin.md)
    with stepwise selection;
    [`reg_predict()`](https://zhjx19.github.io/mathmodels/reference/reg_predict.md),
    [`reg_diagnostics()`](https://zhjx19.github.io/mathmodels/reference/reg_diagnostics.md),
    [`plot_reg_predict()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_predict.md),
    [`plot_reg_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_residuals.md)),
  - lightweight **time series toolkit** built around `ts_df`:
    [`as_ts_df()`](https://zhjx19.github.io/mathmodels/reference/as_ts_df.md)
    for conversion,
    [`complete_ts_df()`](https://zhjx19.github.io/mathmodels/reference/complete_ts_df.md)
    /
    [`impute_ts_df()`](https://zhjx19.github.io/mathmodels/reference/impute_ts_df.md)
    for gap handling,
    [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
    /
    [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md)
    for Box-Cox / log / differencing workflows,
    [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md)
    for stationarity tests,
    [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md)
    for STL decomposition,
    [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
    [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
    and
    [`ts_arimax()`](https://zhjx19.github.io/mathmodels/reference/ts_arimax.md)
    for modeling,
    [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md)
    for unified forecasting, plus
    [`plot_ts()`](https://zhjx19.github.io/mathmodels/reference/plot_ts.md),
    [`plot_ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_forecast.md),
    [`plot_ts_acf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_acf.md),
    [`plot_ts_pacf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_pacf.md),
    [`plot_ts_stl()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_stl.md),
    and
    [`plot_ts_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_residuals.md).
- **interpolation and curve fitting** — piecewise linear, cubic spline,
  polynomial, and Hermite interpolation
  ([`interp_linear()`](https://zhjx19.github.io/mathmodels/reference/interp_linear.md),
  [`interp_spline()`](https://zhjx19.github.io/mathmodels/reference/interp_spline.md),
  [`interp_poly()`](https://zhjx19.github.io/mathmodels/reference/interp_poly.md),
  [`interp_hermite()`](https://zhjx19.github.io/mathmodels/reference/interp_hermite.md)),
  plus curve fitting
  ([`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md),
  [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md),
  [`growth_fit()`](https://zhjx19.github.io/mathmodels/reference/growth_fit.md))
  with automatic starting values for logistic, Gompertz, and
  Michaelis-Menten nonlinear models.
- **Differential Equation Models** — String-formula
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  for arbitrary ODE systems; ready-to-use population models (Malthus,
  Logistic), epidemic compartment models (SI, SIS, SIR, SEIR), and
  Lotka–Volterra predator–prey model, all with a unified `init` +
  `params` interface. Epidemic visualization
  ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
  [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
  and metrics
  ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md)
  for R0, peak, attack rate).
- **Tidyverse Integration** — Seamlessly works with `|>` and
  `dplyr`/`ggplot2` for smooth data manipulation and batch processing.

## Installation

Once released on CRAN, install `mathmodels` with:

``` r

install.packages("mathmodels")
```

You can install the latest development version directly from
[GitHub](https://github.com/zhjx19/mathmodels):

``` r

remotes::install_github("zhjx19/mathmodels")
```

## Getting Started

``` r

library(mathmodels)

# --- AHP (Analytic Hierarchy Process) ---
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
AHP(A)

# --- Epidemic compartment modeling ---
result = ode_sir(
  init   = c(S = 990, I = 10, R = 0),
  params = c(beta = 0.15, gamma = 0.1),  # R0 = 1.5
  times  = seq(0, 100, by = 0.5))

# Visualize
plot_compartments(result, compartments = c("S", "I", "R"))

# Compute epidemic metrics
metrics = epi_metrics(result, beta = 0.15, gamma = 0.1, N = 1000)
metrics$R0         # basic reproduction number

# --- Time series ---
x = as_ts_df(log(AirPassengers))
fit = ts_sarima(x)
fc  = ts_forecast(fit, h = 12)
plot_ts_forecast(x, fc)           # forecast with CIs

# --- Statistical inference (tidy in, tidy out) ---
res = stat_t_test(sleep, .cols = extra, group = group)
res                               # one tidy row: statistic, p, CI, Cohen's d
export_table(res, output = "latex", caption = "Two-sample t-test")

# --- Multivariate statistics ---
r = mv_pca(mtcars, cyl:carb)
plot_mv_scree(r, type = "cum")    # cumulative variance scree
```

## Learning More

For detailed documentation, tutorials, and in-depth examples on using
the `mathmodels` package, please refer to the comprehensive online
manual:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling
(Online Book, In Chinese)](https://zhjx19.github.io/mathmodels-book/)**

This online book is the definitive guide to the package’s
functionalities. Currently implemented modules include:

- **Differential equation models**: Malthus, Logistic, SI, SIS, SIR,
  SEIR, Lotka–Volterra with
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  and `ode_*()` functions; epidemic visualization
  ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
  [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
  and metrics
  ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md))
- **Statistical inference** (0.0.13): tidy tests —
  [`stat_describe()`](https://zhjx19.github.io/mathmodels/reference/stat_describe.md),
  [`stat_t_test()`](https://zhjx19.github.io/mathmodels/reference/stat_t_test.md),
  [`stat_wilcox_test()`](https://zhjx19.github.io/mathmodels/reference/stat_wilcox_test.md),
  [`stat_anova()`](https://zhjx19.github.io/mathmodels/reference/stat_anova.md),
  [`stat_kruskal()`](https://zhjx19.github.io/mathmodels/reference/stat_kruskal.md),
  [`stat_ancova()`](https://zhjx19.github.io/mathmodels/reference/stat_ancova.md),
  [`stat_cor()`](https://zhjx19.github.io/mathmodels/reference/stat_cor.md),
  [`stat_chisq()`](https://zhjx19.github.io/mathmodels/reference/stat_chisq.md),
  [`stat_normality()`](https://zhjx19.github.io/mathmodels/reference/stat_normality.md),
  [`stat_hov_test()`](https://zhjx19.github.io/mathmodels/reference/stat_hov_test.md)
  — all with `.by` slice grouping and effect sizes; generic frameworks
  [`stat_bootstrap()`](https://zhjx19.github.io/mathmodels/reference/stat_bootstrap.md),
  [`stat_perm_test()`](https://zhjx19.github.io/mathmodels/reference/stat_perm_test.md);
  [`export_table()`](https://zhjx19.github.io/mathmodels/reference/export_table.md)
  for Word/LaTeX three-line tables
- **Multivariate statistics** (0.0.12):
  [`mv_pca()`](https://zhjx19.github.io/mathmodels/reference/mv_pca.md),
  [`mv_efa()`](https://zhjx19.github.io/mathmodels/reference/mv_efa.md),
  [`mv_hclust()`](https://zhjx19.github.io/mathmodels/reference/mv_hclust.md),
  [`mv_kmeans()`](https://zhjx19.github.io/mathmodels/reference/mv_kmeans.md),
  [`mv_cluster_metrics()`](https://zhjx19.github.io/mathmodels/reference/mv_cluster_metrics.md),
  [`mv_discrim()`](https://zhjx19.github.io/mathmodels/reference/mv_discrim.md),
  [`mv_corresp()`](https://zhjx19.github.io/mathmodels/reference/mv_corresp.md),
  [`mv_mds()`](https://zhjx19.github.io/mathmodels/reference/mv_mds.md),
  [`mv_cancor()`](https://zhjx19.github.io/mathmodels/reference/mv_cancor.md),
  [`mv_pcor()`](https://zhjx19.github.io/mathmodels/reference/mv_pcor.md)
  with `plot_mv_*()` visuals
- **Data cleaning** (0.0.14):
  [`na_summary()`](https://zhjx19.github.io/mathmodels/reference/na_summary.md),
  [`na_mcar_test()`](https://zhjx19.github.io/mathmodels/reference/na_mcar_test.md)
  (Little’s MCAR test),
  [`plot_na_heatmap()`](https://zhjx19.github.io/mathmodels/reference/plot_na_heatmap.md)
  /
  [`plot_na_bar()`](https://zhjx19.github.io/mathmodels/reference/plot_na_bar.md);
  imputation via
  [`impute()`](https://zhjx19.github.io/mathmodels/reference/impute.md)
  (group mean/median/mode, linear/spline, constant), the unified
  model-based framework
  [`impute_model()`](https://zhjx19.github.io/mathmodels/reference/impute_model.md)
  (kNN, linear regression, decision trees), and multiple imputation
  [`impute_multiple()`](https://zhjx19.github.io/mathmodels/reference/impute_multiple.md)
  (FCS chained equations) with Rubin’s-rules pooling
  [`rubin_pool()`](https://zhjx19.github.io/mathmodels/reference/rubin_pool.md)
  /
  [`pool_fit()`](https://zhjx19.github.io/mathmodels/reference/pool_fit.md);
  outlier handling with
  [`outlier()`](https://zhjx19.github.io/mathmodels/reference/outlier.md)
  (detect / set NA / cap)
- **Pivot table** (0.0.14):
  [`pivot_table()`](https://zhjx19.github.io/mathmodels/reference/pivot_table.md)
  for Excel-style aggregation and spreading
- **Time series**:
  [`ts_df()`](https://zhjx19.github.io/mathmodels/reference/ts_df.md),
  [`as_ts_df()`](https://zhjx19.github.io/mathmodels/reference/as_ts_df.md),
  [`validate_ts_df()`](https://zhjx19.github.io/mathmodels/reference/validate_ts_df.md),
  [`complete_ts_df()`](https://zhjx19.github.io/mathmodels/reference/complete_ts_df.md),
  [`impute_ts_df()`](https://zhjx19.github.io/mathmodels/reference/impute_ts_df.md),
  [`drop_na_ts_df()`](https://zhjx19.github.io/mathmodels/reference/drop_na_ts_df.md),
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md),
  [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md),
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
  [`ts_arimax()`](https://zhjx19.github.io/mathmodels/reference/ts_arimax.md),
  [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md),
  [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md),
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md);
  visualization:
  [`plot_ts()`](https://zhjx19.github.io/mathmodels/reference/plot_ts.md),
  [`plot_ts_acf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_acf.md),
  [`plot_ts_pacf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_pacf.md),
  [`plot_ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_forecast.md),
  [`plot_ts_stl()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_stl.md),
  [`plot_ts_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_residuals.md)
- **Regression prediction**:
  [`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
  [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md),
  [`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
  [`reg_negbin()`](https://zhjx19.github.io/mathmodels/reference/reg_negbin.md)
  with stepwise selection;
  [`reg_predict()`](https://zhjx19.github.io/mathmodels/reference/reg_predict.md),
  [`reg_diagnostics()`](https://zhjx19.github.io/mathmodels/reference/reg_diagnostics.md),
  [`plot_reg_predict()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_predict.md),
  [`plot_reg_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_reg_residuals.md)
- **Interpolation & curve fitting**:
  [`interp_linear()`](https://zhjx19.github.io/mathmodels/reference/interp_linear.md),
  [`interp_spline()`](https://zhjx19.github.io/mathmodels/reference/interp_spline.md),
  [`interp_poly()`](https://zhjx19.github.io/mathmodels/reference/interp_poly.md),
  [`interp_hermite()`](https://zhjx19.github.io/mathmodels/reference/interp_hermite.md)
  for interpolation;
  [`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md),
  [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md),
  [`growth_fit()`](https://zhjx19.github.io/mathmodels/reference/growth_fit.md)
  for curve fitting
- **Markov chain prediction**:
  [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md)
  and
  [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)
- Indicator data preprocessing
- AHP, Entropy weighting, CRITIC, CV, PCA and factor analysis weighting
  (`weight_*()` family,
  incl. [`weight_fa()`](https://zhjx19.github.io/mathmodels/reference/weight_fa.md))
- Weight combination techniques
  ([`combine_weights()`](https://zhjx19.github.io/mathmodels/reference/combine_weights.md))
- TOPSIS, Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE) with the
  `mf_*` membership-function family
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree
