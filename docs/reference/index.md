# Package index

## Evaluation Algorithms

Evaluation, weighting, fuzzy evaluation, DEA, inequality, regional
economics, and system evaluation methods.

## Data Preprocessing

Scaling, normalization, and indicator transformation utilities.

- [`standardize()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`normalize()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_middle()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_interval()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_extreme()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_initial()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_mean()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`to_positive()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  : Preprocessing Functions for Data Normalization and Standardization

## Weighting Methods

Subjective, objective, and combined weighting methods for
multi-indicator evaluation.

- [`AHP()`](https://zhjx19.github.io/mathmodels/reference/AHP.md) : AHP:
  Analytic Hierarchy Process
- [`combine_weights()`](https://zhjx19.github.io/mathmodels/reference/combine_weights.md)
  : Combine Subjective and Objective Weights
- [`critic_weight()`](https://zhjx19.github.io/mathmodels/reference/critic_weight.md)
  : CRITIC Weight Method
- [`cv_weight()`](https://zhjx19.github.io/mathmodels/reference/cv_weight.md)
  : Coefficient of Variation Weighting
- [`entropy_weight()`](https://zhjx19.github.io/mathmodels/reference/entropy_weight.md)
  : Entropy Weight Method
- [`linear_sum()`](https://zhjx19.github.io/mathmodels/reference/linear_sum.md)
  : Linear weighted synthesis
- [`pca_weight()`](https://zhjx19.github.io/mathmodels/reference/pca_weight.md)
  : PCA-Based Weighting Method

## Evaluation Models

Multi-criteria ranking and grey relational evaluation methods.

- [`topsis()`](https://zhjx19.github.io/mathmodels/reference/topsis.md)
  : TOPSIS Method for Multi-Criteria Decision Making
- [`rank_sum_ratio()`](https://zhjx19.github.io/mathmodels/reference/rank_sum_ratio.md)
  : Rank Sum Ratio (RSR) Evaluation
- [`grey_corr()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  [`grey_corr_topsis()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  : Grey Relational Analysis Functions

## Fuzzy Evaluation

Membership functions, fuzzy comprehensive evaluation, and
defuzzification.

- [`tri_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`trap_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gauss_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gbell_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gauss2mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`sigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`dsigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`psigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`z_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`pi_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`s_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`plot_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  : Membership Functions for Fuzzy Logic
- [`compute_mf_funs()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md)
  [`compute_mf()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md)
  : Compute fuzzy membership vector and return corresponding membership
  functions.
- [`fuzzy_eval()`](https://zhjx19.github.io/mathmodels/reference/fuzzy_eval.md)
  : Fuzzy Comprehensive Evaluation
- [`defuzzify()`](https://zhjx19.github.io/mathmodels/reference/defuzzify.md)
  : Defuzzification Methods for Fuzzy Comprehensive Evaluation

## Data Envelopment Analysis

DEA, SBM, super-efficiency, and Malmquist productivity models.

- [`basic_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`super_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`basic_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`super_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`malmquist()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  : DEA efficiency analysis

## Inequality Measures

Gini and Theil indices for individual, grouped, and nested data.

- [`gini0()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`gini()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil0()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil0_g()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g2_cross()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g2_nest()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  : Inequality Indices

## Regional Economics

Regional and industrial concentration indices.

- [`LQ()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  [`HHI()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  [`EG()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  : Regional Economics Functions

## System Evaluation

Coupling coordination and obstacle degree analysis for multi-indicator
systems.

- [`coupling_degree()`](https://zhjx19.github.io/mathmodels/reference/system_evaluation.md)
  [`obstacle_degree()`](https://zhjx19.github.io/mathmodels/reference/system_evaluation.md)
  : System Evaluation Functions for Coupling and Obstacle Analysis

## Prediction Algorithms

Grey prediction, Markov prediction, and time series models.

## Grey Prediction Models

Grey forecasting models and utilities for combining predictions.

- [`GM11()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`GM1N()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`DGM21()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`verhulst()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  : Grey Prediction Models
- [`combine_preds()`](https://zhjx19.github.io/mathmodels/reference/combine_preds.md)
  : Combine Multiple Prediction Results

## Markov Prediction Models

Markov chain prediction and Grey-Markov forecasting models.

- [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md)
  : Markov Chain Prediction
- [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)
  : Grey-Markov Prediction Model

## Time Series Models

Time series transformation, diagnostics, modeling, forecasting, and
visualization.

- [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
  : Transform a Time Series for Stationarity
- [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md)
  : Back-Transform Forecasts to the Original Scale
- [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md)
  : Stationarity Tests for a Time Series
- [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md)
  : STL (Seasonal + Trend + Loess) Decomposition
- [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md)
  : ETS (Error, Trend, Seasonality) Exponential Smoothing
- [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md)
  : SARIMA Model Fitting
- [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md)
  : GARCH Variance Modeling
- [`ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima_garch.md)
  : Two-Stage SARIMA-GARCH Joint Model
- [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md)
  : Generate Forecasts from a Fitted Time Series Model
- [`plot_ts()`](https://zhjx19.github.io/mathmodels/reference/plot_ts.md)
  : Plot a Time Series
- [`plot_ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_forecast.md)
  : Plot Historical Series + Forecast with Confidence Bands
- [`plot_ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_sarima_garch.md)
  : Dual-Axis Plot for SARIMA-GARCH: Mean + Conditional Volatility
- [`plot_ts_garch()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_garch.md)
  : GARCH Volatility Plot
- [`plot_ts_residuals()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_residuals.md)
  : Residual Diagnostic Plots
- [`plot_ts_stl()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_stl.md)
  : Plot STL Decomposition Components
- [`plot_ts_acf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_acf.md)
  : ACF and PACF Plots
- [`plot_ts_pacf()`](https://zhjx19.github.io/mathmodels/reference/plot_ts_pacf.md)
  : PACF Plot (convenience alias)

## Differential Equation Algorithms

ODE solvers, dynamic models, epidemic models, metrics, and
visualizations.

## Differential Equation Models

ODE solving, population dynamics, epidemic models, metrics, and
visualizations.

- [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  : General ODE Solver
- [`model_malthus()`](https://zhjx19.github.io/mathmodels/reference/model_malthus.md)
  : Malthusian (Exponential) Growth Model
- [`model_logistic()`](https://zhjx19.github.io/mathmodels/reference/model_logistic.md)
  : Logistic Population Growth Model
- [`model_si()`](https://zhjx19.github.io/mathmodels/reference/model_si.md)
  : SI Epidemic Model
- [`model_sis()`](https://zhjx19.github.io/mathmodels/reference/model_sis.md)
  : SIS Epidemic Model
- [`model_sir()`](https://zhjx19.github.io/mathmodels/reference/model_sir.md)
  : SIR Epidemic Model
- [`model_seir()`](https://zhjx19.github.io/mathmodels/reference/model_seir.md)
  : SEIR Epidemic Model
- [`model_lv()`](https://zhjx19.github.io/mathmodels/reference/model_lv.md)
  : Lotka-Volterra Predator-Prey Model
- [`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md)
  : Extract key epidemic metrics from ODE simulation output
- [`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md)
  : Plot compartment trajectories
- [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md)
  : Plot daily new infections (dI)
- [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md)
  : Phase plot S vs I
- [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md)
  : Plot effective reproduction number R_t

## Data and Utilities

Built-in datasets and data import utilities.

## Datasets

- [`water_quality`](https://zhjx19.github.io/mathmodels/reference/water_quality.md)
  : Water Quality Dataset

## Data I/O Utilities

- [`read_nbs()`](https://zhjx19.github.io/mathmodels/reference/read_nbs.md)
  : Read and Combine National Bureau of Statistics XLS Files
