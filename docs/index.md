# mathmodels: A Powerful R Toolkit for Mathematical Modeling

An R package providing a versatile toolkit for mathematical modeling,
developed as a companion to the book *Mathematical Modeling: Algorithms
and Programming Implementation* (China Machine Press). It focuses on
implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.10)** introduces a regression prediction toolkit
(`reg_lm`, `reg_logistic`, `reg_poisson`, `reg_negbin`) with stepwise
selection, diagnostics, and visualization — on top of the time-series
workflow, epidemic modeling, differential equation models, grey
prediction, Markov chain models, and a rich suite of evaluation methods
(AHP, Entropy, CRITIC, PCA, TOPSIS, Fuzzy, RSR, DEA).

## Key Features

- **Evaluation Models** — AHP, Entropy weighting, CRITIC, PCA weighting,
  TOPSIS, Grey Relational Analysis (GRA), Rank Sum Ratio (RSR), Fuzzy
  Comprehensive Evaluation (FCE), Data Envelopment Analysis
  (CCR/BCC/SBM, Malmquist), plus inequality measures (Gini, Theil
  Index), coupling coordination degree, and obstacle degree models.
- **Prediction Models** — Grey prediction (GM(1,1), GM(1,N), Verhulst),
  Markov chain prediction
  ([`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md),
  [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)),
  and a lightweight **time series toolkit** built around `ts_df`:
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
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md)
  and
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md)
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
result = model_sir(
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
  and `model_*()` functions; epidemic visualization
  ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
  [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
  and metrics
  ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md))
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
- **Markov chain prediction**:
  [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md)
  and
  [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)
- Indicator data preprocessing
- AHP, Entropy weighting, CRITIC, PCA weighting
- Weight combination techniques
- TOPSIS, Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE)
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree
