
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mathmodels:  A Powerful R Toolkit for Mathematical Modeling

<!-- badges: start -->
<!-- badges: end -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mathmodels)](https://cran.r-project.org/package=mathmodels)
[![Online Manual](https://img.shields.io/badge/Online-Manual-blue)](https://zhjx19.github.io/mathmodels-book/)

An R package providing a versatile toolkit for mathematical modeling, developed as a companion to the book *Mathematical Modeling: Algorithms and Programming Implementation* (China Machine Press). It focuses on implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.10)** introduces a regression prediction toolkit (`reg_lm`,
`reg_logistic`, `reg_poisson`, `reg_negbin`) with stepwise selection, diagnostics,
and visualization — on top of the time-series workflow, epidemic modeling,
differential equation models, grey prediction, Markov chain models, and a rich
suite of evaluation methods (AHP, Entropy, CRITIC, PCA, TOPSIS, Fuzzy, RSR, DEA).

## Key Features

*   **Evaluation Models** — AHP, Entropy weighting, CRITIC, PCA weighting, TOPSIS, Grey Relational Analysis (GRA),
    Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE), Data Envelopment Analysis (CCR/BCC/SBM, Malmquist),
    plus inequality measures (Gini, Theil Index), coupling coordination degree, and obstacle degree models.
*   **Prediction Models** — 
    - **Grey prediction** (GM(1,1), GM(1,N), Verhulst), Markov chain prediction (`markov_chain()`,`GM11_markov()`), 
    - **regression prediction** toolkit (`reg_lm()`, `reg_logistic()`, `reg_poisson()`,
      `reg_negbin()` with stepwise selection; `reg_predict()`, `reg_diagnostics()`, `plot_reg_predict()`, `plot_reg_residuals()`), 
    - lightweight **time series toolkit** built around `ts_df`: `as_ts_df()` for conversion,
      `complete_ts_df()` / `impute_ts_df()` for gap handling, `ts_transform()` / `ts_back_transform()` for
      Box-Cox / log / differencing workflows, `ts_test()` for stationarity tests, `ts_stl()` for STL decomposition,
      `ts_ets()` and `ts_sarima()` for modeling, `ts_forecast()` for unified forecasting, plus `plot_ts()`,
      `plot_ts_forecast()`, `plot_ts_acf()`, `plot_ts_pacf()`, `plot_ts_stl()`, and `plot_ts_residuals()`.
*   **interpolation and curve fitting** — piecewise linear, cubic spline, polynomial, and Hermite
    interpolation (`interp_linear()`, `interp_spline()`, `interp_poly()`, `interp_hermite()`), plus curve fitting
    (`poly_fit()`, `curve_fit()`, `growth_fit()`) with automatic starting values for logistic, Gompertz,
    and Michaelis-Menten nonlinear models.
*   **Differential Equation Models** — String-formula `ode_solver()` for arbitrary ODE systems; ready-to-use
    population models (Malthus, Logistic), epidemic compartment models (SI, SIS, SIR, SEIR), and
    Lotka–Volterra predator–prey model, all with a unified `init` + `params` interface.
    Epidemic visualization (`plot_compartments()`, `plot_incidence()`, `plot_phase_si()`,
    `plot_Rt_estimate()`) and metrics (`epi_metrics()` for R0, peak, attack rate).
*   **Tidyverse Integration** — Seamlessly works with `|>` and `dplyr`/`ggplot2` for
    smooth data manipulation and batch processing.

## Installation

Once released on CRAN, install `mathmodels` with:

``` r
install.packages("mathmodels")
```

You can install the latest development version directly from [GitHub](https://github.com/zhjx19/mathmodels):

``` r
remotes::install_github("zhjx19/mathmodels")
```

## Getting Started

```r
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

For detailed documentation, tutorials, and in-depth examples on using the `mathmodels` package, please refer to the comprehensive online manual:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling (Online Book, In Chinese)](https://zhjx19.github.io/mathmodels-book/)**

This online book is the definitive guide to the package's functionalities. Currently implemented modules include:

- **Differential equation models**: Malthus, Logistic, SI, SIS, SIR, SEIR, Lotka–Volterra with `ode_solver()` and `model_*()` functions; epidemic visualization (`plot_compartments()`, `plot_incidence()`, `plot_phase_si()`, `plot_Rt_estimate()`) and metrics (`epi_metrics()`)
- **Time series**: `ts_df()`, `as_ts_df()`, `validate_ts_df()`, `complete_ts_df()`, `impute_ts_df()`, `drop_na_ts_df()`, `ts_transform()`, `ts_back_transform()`, `ts_ets()`, `ts_sarima()`, `ts_stl()`, `ts_test()`, `ts_forecast()`; visualization: `plot_ts()`, `plot_ts_acf()`, `plot_ts_pacf()`, `plot_ts_forecast()`, `plot_ts_stl()`, `plot_ts_residuals()`
- **Regression prediction**: `reg_lm()`, `reg_logistic()`, `reg_poisson()`, `reg_negbin()` with stepwise selection; `reg_predict()`, `reg_diagnostics()`, `plot_reg_predict()`, `plot_reg_residuals()`
- **Interpolation & curve fitting**: `interp_linear()`, `interp_spline()`, `interp_poly()`, `interp_hermite()` for interpolation; `poly_fit()`, `curve_fit()`, `growth_fit()` for curve fitting
- **Markov chain prediction**: `markov_chain()` and `GM11_markov()`
- Indicator data preprocessing
- AHP, Entropy weighting, CRITIC, PCA weighting
- Weight combination techniques
- TOPSIS, Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE)
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree







